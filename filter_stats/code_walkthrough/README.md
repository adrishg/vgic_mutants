# code_tutor: step-through walkthroughs of the pipeline's code

Five interactive pages, one per pipeline file, in the order the pipeline runs them:

| Part | Page | File | What it does |
|---|---|---|---|
| 1 | `pages/config.html` | `config.py` | reads and checks the dataset config |
| 2 | `pages/build_dv.html` | `build_dv.py` | step 1: archives → numbered, superposed models → DVs |
| 3 | `pages/screen.html` | `screen.py` | the topology screen: removes models with the wrong fold |
| 4 | `pages/region_dv.html` | `region_dv.py` | step 2: the DVs of each region, on the kept models |
| 5 | `pages/stats.html` | `stats.py` | step 3: mixed model and comparisons with the control |

Each page ran the real file, unchanged, on one small made-up dataset, and shows every statement with the values
it produced. Step into a function to follow it line by line, or step over it to see only its result. Where a
table's rows, columns or positions stand for something the numbers don't show (which model, chain, residue or
comparison), a second box labels them. Each step's output files are the next step's input.

## Opening the pages

GitHub shows `.html` files as source code, so open them in a browser instead:

1. Download this folder (or clone the repository).
2. Open `index.html` in a browser and pick a walkthrough.

Keep `index.html` and `pages/` together so the links between them work. The pages work offline; only the fonts
come from the web (without a connection a system font is used).

## Notes

- The pages show the pipeline's code as it was when they were built. If the code has changed since, line
  numbers and values may no longer match the current files.
- The toy dataset is made up: a 43-residue, four-chain channel with two conditions (WT and R18Q), two networks
  and three seeds. It is small enough to read on screen, not realistic, and contains nothing from real datasets.
