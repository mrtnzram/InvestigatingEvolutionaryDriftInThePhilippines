"""Build and open a standalone HTML reviewer for Philippine loanword labels.

Run from the project root with:

    python python/cognate/loanword_reviewer.py

The reviewer reads output/PH_df.csv and output/PH_loan.csv, preserves review
progress in the browser, and downloads the confirmed selections as a CSV.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import webbrowser
from collections import OrderedDict
from pathlib import Path

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "output"
DEFAULT_PH_DF = DEFAULT_OUTPUT_DIR / "PH_df.csv"
DEFAULT_PH_LOAN = DEFAULT_OUTPUT_DIR / "PH_loan.csv"
DEFAULT_HTML = DEFAULT_OUTPUT_DIR / "loanword_reviewer.html"

METADATA_COLUMNS = {
    "language_id",
    "language",
    "latitude",
    "longitude",
    "glottocode",
    "author",
    "number_of_entries",
    "source_url",
}


def split_cell(value: object, separator: str) -> list[str]:
    if pd.isna(value):
        return []
    return [part.strip() for part in str(value).split(separator) if part.strip()]


def append_unique(values: list[str], value: str) -> None:
    key = value.casefold()
    if all(existing.casefold() != key for existing in values):
        values.append(value)


def build_concepts(
    ph_df: pd.DataFrame,
    ph_loan: pd.DataFrame,
    word_columns: list[str],
) -> list[dict]:
    missing = [column for column in word_columns if column not in ph_loan.columns]
    if missing:
        raise ValueError(f"PH_loan is missing {len(missing)} word columns: {missing[:5]}")

    if "language_id" not in ph_df.columns or "language_id" not in ph_loan.columns:
        raise ValueError("Both inputs must contain a language_id column.")

    loan_lookup: dict[str, pd.Series] = {}
    for _, row in ph_loan.iterrows():
        language_key = str(row["language_id"])
        if language_key in loan_lookup:
            raise ValueError(f"PH_loan has a duplicate language_id: {language_key}")
        loan_lookup[language_key] = row

    concepts: list[dict] = []

    for english_word in word_columns:
        unique_forms: list[str] = []
        evidence_by_form: OrderedDict[str, dict] = OrderedDict()

        for _, row in ph_df.iterrows():
            forms = split_cell(row.get(english_word), ";")
            if not forms:
                continue

            for form in forms:
                append_unique(unique_forms, form)

            language_key = str(row["language_id"])
            loan_row = loan_lookup.get(language_key)
            annotations = (
                split_cell(loan_row.get(english_word), ",")
                if loan_row is not None
                else []
            )
            if not annotations:
                continue

            language_name = str(row.get("language", "")).strip()
            for form in forms:
                evidence = evidence_by_form.setdefault(
                    form,
                    {"word": form, "annotations": [], "languages": []},
                )
                for annotation in annotations:
                    append_unique(evidence["annotations"], annotation)
                if language_name and language_name.lower() != "nan":
                    append_unique(evidence["languages"], language_name)

        concepts.append(
            {
                "english_word": english_word,
                "n_unique": len(unique_forms),
                "unique_words": unique_forms,
                "evidence": list(evidence_by_form.values()),
                "initial_selected": list(evidence_by_form.keys()),
            }
        )

    return sorted(
        concepts,
        key=lambda concept: (-concept["n_unique"], concept["english_word"].casefold()),
    )


def build_html(concepts: list[dict]) -> str:
    compact_json = json.dumps(concepts, ensure_ascii=False, separators=(",", ":"))
    compact_json = compact_json.replace("</", "<\\/")
    fingerprint = hashlib.sha256(compact_json.encode("utf-8")).hexdigest()[:16]

    template = r'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Philippine Loanword Reviewer</title>
  <style>
    :root {
      color-scheme: light dark;
      --bg: light-dark(#f4f6f8, #111418);
      --panel: light-dark(#ffffff, #1a1f25);
      --panel-soft: light-dark(#f7f9fb, #20262d);
      --text: light-dark(#17202a, #ecf1f5);
      --muted: light-dark(#617080, #a6b0ba);
      --border: light-dark(#d8dee5, #3a424b);
      --accent: light-dark(#1f6feb, #63a4ff);
      --accent-soft: light-dark(#eaf2ff, #183152);
      --success: light-dark(#16803c, #52c878);
      --success-soft: light-dark(#eaf7ee, #173b25);
      --shadow: light-dark(0 14px 34px rgba(31, 43, 55, .10), 0 14px 34px rgba(0, 0, 0, .28));
      font-family: Inter, ui-sans-serif, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }

    * { box-sizing: border-box; }
    body { margin: 0; background: var(--bg); color: var(--text); }
    button, input { font: inherit; }

    .app {
      width: min(1180px, calc(100% - 32px));
      margin: 24px auto;
    }

    .topbar {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 20px;
      margin-bottom: 14px;
    }

    .progress-copy { color: var(--muted); font-size: 14px; }
    .progress-track {
      flex: 1;
      height: 8px;
      overflow: hidden;
      border-radius: 999px;
      background: var(--border);
    }
    .progress-fill {
      width: 0;
      height: 100%;
      border-radius: inherit;
      background: var(--success);
      transition: width .2s ease;
    }

    .workspace {
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 18px;
      box-shadow: var(--shadow);
      overflow: hidden;
    }

    .concept-header {
      display: flex;
      justify-content: space-between;
      align-items: center;
      gap: 18px;
      padding: 22px 24px;
      border-bottom: 1px solid var(--border);
    }
    .concept-header h1 { margin: 0; font-size: clamp(24px, 4vw, 34px); font-weight: 500; }
    .count-badge {
      flex: none;
      padding: 7px 11px;
      border-radius: 999px;
      color: var(--accent);
      background: var(--accent-soft);
      font-size: 14px;
      font-weight: 500;
    }

    .columns {
      display: grid;
      grid-template-columns: minmax(0, 1.35fr) minmax(300px, .8fr);
      min-height: 430px;
    }
    .pane { padding: 22px 24px 26px; }
    .pane + .pane { border-left: 1px solid var(--border); background: var(--panel-soft); }
    .pane-title { margin: 0 0 4px; font-size: 17px; font-weight: 500; }
    .pane-subtitle { margin: 0 0 18px; color: var(--muted); font-size: 14px; }

    .word-grid {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
      gap: 10px;
    }
    .word-option {
      display: flex;
      align-items: center;
      gap: 10px;
      min-height: 44px;
      padding: 9px 11px;
      border: 1px solid var(--border);
      border-radius: 10px;
      background: var(--panel);
      cursor: pointer;
      overflow-wrap: anywhere;
    }
    .word-option:has(input:checked) {
      border-color: var(--accent);
      background: var(--accent-soft);
    }
    .word-option input { width: 17px; height: 17px; accent-color: var(--accent); flex: none; }

    .evidence-list { display: grid; gap: 10px; }
    .evidence-item {
      padding: 12px;
      border: 1px solid var(--border);
      border-radius: 10px;
      background: var(--panel);
    }
    .evidence-item.selected { border-color: var(--success); background: var(--success-soft); }
    .evidence-word { font-weight: 500; overflow-wrap: anywhere; }
    .annotations { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px; }
    .annotation {
      padding: 3px 7px;
      border-radius: 999px;
      background: var(--accent-soft);
      color: var(--accent);
      font-size: 12px;
    }
    .languages { margin-top: 8px; color: var(--muted); font-size: 12px; overflow-wrap: anywhere; }
    .empty { color: var(--muted); font-size: 14px; }

    .footer {
      display: grid;
      grid-template-columns: auto 1fr auto;
      align-items: center;
      gap: 12px;
      padding: 16px 20px;
      border-top: 1px solid var(--border);
    }
    .footer-center { display: flex; justify-content: center; gap: 10px; }
    .button {
      min-height: 40px;
      padding: 8px 14px;
      border: 1px solid var(--border);
      border-radius: 9px;
      background: var(--panel);
      color: var(--text);
      cursor: pointer;
    }
    .button:hover:not(:disabled) { border-color: var(--accent); }
    .button.primary { border-color: var(--accent); background: var(--accent); color: white; }
    .button.export { border-color: var(--success); color: var(--success); }
    .button:disabled { opacity: .45; cursor: not-allowed; }

    .saved-note { margin-top: 10px; text-align: right; color: var(--muted); font-size: 12px; }

    @media (max-width: 800px) {
      .app { width: min(100% - 20px, 680px); margin: 10px auto; }
      .topbar { align-items: flex-start; flex-direction: column; gap: 8px; }
      .progress-track { width: 100%; }
      .columns { grid-template-columns: 1fr; }
      .pane + .pane { border-left: 0; border-top: 1px solid var(--border); }
      .footer { grid-template-columns: 1fr 1fr; }
      .footer-center { grid-column: 1 / -1; grid-row: 1; }
      .footer > .button:first-child { grid-column: 1; }
      .footer > .button:last-child { grid-column: 2; }
    }
  </style>
</head>
<body>
  <main class="app">
    <div class="topbar">
      <div class="progress-copy" id="progressCopy"></div>
      <div class="progress-track" role="progressbar" aria-label="Review progress" aria-valuemin="0" aria-valuemax="100" aria-valuenow="0" id="progressTrack">
        <div class="progress-fill" id="progressFill"></div>
      </div>
    </div>

    <section class="workspace">
      <header class="concept-header">
        <h1 id="conceptTitle"></h1>
        <div class="count-badge" id="uniqueCount"></div>
      </header>

      <div class="columns">
        <section class="pane" aria-labelledby="candidateHeading">
          <h2 class="pane-title" id="candidateHeading">Candidate words</h2>
          <p class="pane-subtitle">Select every form that should be treated as a loan.</p>
          <div class="word-grid" id="wordGrid"></div>
        </section>

        <aside class="pane" aria-labelledby="evidenceHeading">
          <h2 class="pane-title" id="evidenceHeading">Marked evidence</h2>
          <p class="pane-subtitle">Forms initially selected from the annotation table.</p>
          <div class="evidence-list" id="evidenceList"></div>
        </aside>
      </div>

      <footer class="footer">
        <button class="button" id="previousButton" type="button">Previous</button>
        <div class="footer-center">
          <button class="button primary" id="confirmButton" type="button">Confirm &amp; next</button>
        </div>
        <button class="button" id="nextButton" type="button">Next</button>
      </footer>
    </section>

    <div class="saved-note" id="savedNote">Selections save automatically in this browser.</div>
    <div style="display:flex;justify-content:flex-end;margin-top:10px">
      <button class="button export" id="exportButton" type="button" disabled>Download confirmed CSV</button>
    </div>
  </main>

  <script>
    const concepts = __CONCEPTS__;
    const storageKey = "philippine-loanword-review-__FINGERPRINT__";
    const saved = JSON.parse(localStorage.getItem(storageKey) || "null");

    const state = {
      current: 0,
      selections: {},
      confirmed: {},
    };

    for (const concept of concepts) {
      const key = concept.english_word;
      const restored = saved?.selections?.[key];
      state.selections[key] = new Set(Array.isArray(restored) ? restored : concept.initial_selected);
      state.confirmed[key] = Boolean(saved?.confirmed?.[key]);
    }
    if (Number.isInteger(saved?.current)) {
      state.current = Math.max(0, Math.min(saved.current, concepts.length - 1));
    }

    const elements = {
      conceptTitle: document.getElementById("conceptTitle"),
      uniqueCount: document.getElementById("uniqueCount"),
      wordGrid: document.getElementById("wordGrid"),
      evidenceList: document.getElementById("evidenceList"),
      previousButton: document.getElementById("previousButton"),
      nextButton: document.getElementById("nextButton"),
      confirmButton: document.getElementById("confirmButton"),
      exportButton: document.getElementById("exportButton"),
      progressCopy: document.getElementById("progressCopy"),
      progressFill: document.getElementById("progressFill"),
      progressTrack: document.getElementById("progressTrack"),
      savedNote: document.getElementById("savedNote"),
    };

    function persist() {
      const selections = {};
      for (const [key, values] of Object.entries(state.selections)) {
        selections[key] = [...values];
      }
      localStorage.setItem(storageKey, JSON.stringify({
        current: state.current,
        selections,
        confirmed: state.confirmed,
      }));
      elements.savedNote.textContent = "Saved automatically.";
    }

    function render() {
      const concept = concepts[state.current];
      const selected = state.selections[concept.english_word];
      const confirmedCount = Object.values(state.confirmed).filter(Boolean).length;
      const percent = concepts.length ? Math.round(100 * confirmedCount / concepts.length) : 0;

      elements.conceptTitle.textContent = concept.english_word;
      elements.uniqueCount.textContent = `${concept.n_unique} unique`;
      elements.progressCopy.textContent = `Concept ${state.current + 1} of ${concepts.length} · ${confirmedCount} confirmed`;
      elements.progressFill.style.width = `${percent}%`;
      elements.progressTrack.setAttribute("aria-valuenow", String(percent));

      elements.wordGrid.replaceChildren();
      if (!concept.unique_words.length) {
        const empty = document.createElement("div");
        empty.className = "empty";
        empty.textContent = "No recorded forms.";
        elements.wordGrid.append(empty);
      }

      for (const word of concept.unique_words) {
        const label = document.createElement("label");
        label.className = "word-option";
        const input = document.createElement("input");
        input.type = "checkbox";
        input.checked = selected.has(word);
        input.setAttribute("aria-label", `Mark ${word} as a loan`);
        input.addEventListener("change", () => {
          if (input.checked) selected.add(word);
          else selected.delete(word);
          state.confirmed[concept.english_word] = false;
          persist();
          render();
        });
        const text = document.createElement("span");
        text.textContent = word;
        label.append(input, text);
        elements.wordGrid.append(label);
      }

      elements.evidenceList.replaceChildren();
      if (!concept.evidence.length) {
        const empty = document.createElement("div");
        empty.className = "empty";
        empty.textContent = "No annotation evidence was detected for this concept.";
        elements.evidenceList.append(empty);
      }

      for (const evidence of concept.evidence) {
        const item = document.createElement("article");
        item.className = `evidence-item${selected.has(evidence.word) ? " selected" : ""}`;
        const word = document.createElement("div");
        word.className = "evidence-word";
        word.textContent = evidence.word;
        const annotations = document.createElement("div");
        annotations.className = "annotations";
        for (const annotation of evidence.annotations) {
          const chip = document.createElement("span");
          chip.className = "annotation";
          chip.textContent = annotation;
          annotations.append(chip);
        }
        item.append(word, annotations);
        if (evidence.languages.length) {
          const languages = document.createElement("div");
          languages.className = "languages";
          languages.textContent = evidence.languages.join("; ");
          item.append(languages);
        }
        elements.evidenceList.append(item);
      }

      elements.previousButton.disabled = state.current === 0;
      elements.nextButton.disabled = state.current === concepts.length - 1;
      elements.confirmButton.textContent = state.confirmed[concept.english_word]
        ? "Confirmed ✓"
        : (state.current === concepts.length - 1 ? "Confirm" : "Confirm & next");
      elements.exportButton.disabled = confirmedCount !== concepts.length;
    }

    function goTo(index) {
      state.current = Math.max(0, Math.min(index, concepts.length - 1));
      persist();
      render();
      window.scrollTo({ top: 0, behavior: "smooth" });
    }

    elements.previousButton.addEventListener("click", () => goTo(state.current - 1));
    elements.nextButton.addEventListener("click", () => goTo(state.current + 1));
    elements.confirmButton.addEventListener("click", () => {
      const concept = concepts[state.current];
      state.confirmed[concept.english_word] = true;
      persist();
      if (state.current < concepts.length - 1) goTo(state.current + 1);
      else render();
    });

    function csvCell(value) {
      const text = String(value ?? "");
      return `"${text.replaceAll('"', '""')}"`;
    }

    elements.exportButton.addEventListener("click", () => {
      const rows = [["english_word", "detected_loan_words"]];
      for (const concept of concepts) {
        rows.push([
          concept.english_word,
          [...state.selections[concept.english_word]].join("; "),
        ]);
      }
      const csv = rows.map(row => row.map(csvCell).join(",")).join("\n");
      const blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = "confirmed_philippine_loanwords.csv";
      document.body.append(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    });

    render();
  </script>
</body>
</html>
'''

    return template.replace("__CONCEPTS__", compact_json).replace(
        "__FINGERPRINT__", fingerprint
    )


def create_reviewer(
    ph_df_path: Path,
    ph_loan_path: Path,
    output_path: Path,
) -> Path:
    ph_df = pd.read_csv(ph_df_path)
    ph_loan = pd.read_csv(ph_loan_path)
    word_columns = [
        column for column in ph_df.columns if column not in METADATA_COLUMNS
    ]
    if not word_columns:
        raise ValueError("No word columns were found in PH_df.")

    concepts = build_concepts(ph_df, ph_loan, word_columns)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(build_html(concepts), encoding="utf-8")
    return output_path.resolve()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a local interactive reviewer for Philippine loanwords."
    )
    parser.add_argument("--ph-df", type=Path, default=DEFAULT_PH_DF)
    parser.add_argument("--ph-loan", type=Path, default=DEFAULT_PH_LOAN)
    parser.add_argument("--output", type=Path, default=DEFAULT_HTML)
    parser.add_argument(
        "--no-open",
        action="store_true",
        help="Create the HTML without opening it in the default browser.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_path = create_reviewer(args.ph_df, args.ph_loan, args.output)
    print(f"Created {output_path}")
    if not args.no_open:
        webbrowser.open(output_path.as_uri())


if __name__ == "__main__":
    main()
