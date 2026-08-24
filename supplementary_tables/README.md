# Supplementary table CSV exports

These files were exported on 24 August 2026 from the five populated tabs in
`Supplementary tables.gsheet`:

| Google Sheet tab | CSV file | Parsed rows | Columns |
|---|---|---:|---:|
| `meta_data` | `meta_data.csv` | 48 | 4 |
| `Table_S1` | `Table_S1.csv` | 31 | 17 |
| `Table_S2` | `Table_S2.csv` | 168 | 11 |
| `Table_S3` | `Table_S3.csv` | 1,155 | 10 |
| `Table_S4` | `Table_S4.csv` | 675 | 11 |

The four table tabs retain their caption in row 1, a blank row 2 and the table header
in row 3. `meta_data.csv` is header-first. This preserves the visible Google Sheet
content. Header-first machine outputs remain in `revision/results/supplementary/`.

Table S1 contains the evidence-base characteristics and Table S2 contains the reported
metrics. The script-generated files in `revision/results/supplementary/` use the same
numbering.
