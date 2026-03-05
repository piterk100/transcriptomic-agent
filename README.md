# Transcriptomic Agent

Autonomiczny agent AI do odkrywania powiązań transkryptomicznych w wielu datasetach jednocześnie.

## Uruchomienie

```bash
npm install
npm start
```

Otwórz [http://localhost:3000](http://localhost:3000)

## Struktura projektu

```
src/
├── utils/
│   ├── math.js          # funkcje matematyczne (t-test, BH, korelacja, Cohen's d…)
│   └── parser.js        # parsowanie CSV, macierzy ekspresji, metadanych
│
├── tools/
│   ├── singleDataset.js # narzędzia dla jednego datasetu
│   ├── crossDataset.js  # narzędzia cross-dataset (IRM, meta-DE, rewiring…)
│   └── index.js         # rejestr narzędzi + sandbox + summarizeResult
│
├── agent/
│   ├── systemPrompt.js  # budowanie system promptu agenta
│   └── runner.js        # pętla agenta (LLM → wybór narzędzia → wynik → loop)
│
├── components/
│   ├── DatasetSlot.jsx  # komponent slotu do wgrywania plików
│   └── LogEntry.jsx     # komponent wpisu w logu agenta
│
├── App.jsx              # główny komponent, UI
└── index.js             # entry point
```

## Jak dodać nowe narzędzie

1. Napisz funkcję w `src/tools/singleDataset.js` lub `src/tools/crossDataset.js`
   - Przyjmuje `{ datasets, ...params }`
   - Zwraca plain JSON object z wynikami i polem `interpretation`

2. Dodaj wpis do `summarizeResult()` w `src/tools/index.js`

3. Zaktualizuj `buildSystemPrompt()` w `src/agent/systemPrompt.js` — dopisz narzędzie do listy

Agent automatycznie je zobaczy i będzie mógł wywoływać.

## Narzędzia

### Single-dataset
| Narzędzie | Opis |
|---|---|
| `dataset_summary` | Przegląd wszystkich wczytanych datasetów |
| `top_variable_genes` | Top N genów o największej wariancji |
| `differential_expression` | DE między dwoma grupami (t-test + BH) |
| `gene_expression_by_group` | Średnia ekspresja genów w każdej grupie |
| `nonlinear_rule` | Reguła AND: GEN_A > próg AND GEN_B <= próg → grupa |
| `contextual_modules` | Ko-ekspresja zależna od kontekstu (mediacja przez gen) |
| `pathway_enrichment` | Wzbogacenie w 12 pathway (Hallmarks + KEGG) |
| `batch_detection` | Wykrywanie batch artefaktów w odkrytych osiach |
| `subgroup_discovery` | Podgrupy wewnątrz jednej grupy (PCA power iteration) |
| `gene_network_hub` | Huby sieci ko-ekspresji |

### Cross-dataset
| Narzędzie | Opis |
|---|---|
| `cross_dataset_de` | Geny DE w tym samym kierunku we WSZYSTKICH datasetach |
| `cross_dataset_correlation` | Replikacja modułu ko-ekspresji między kohortami |
| `invariant_axis` | Oś biologiczna inwariantna (IRM-style, Cohen's d stability) |
| `cross_dataset_rewiring` | Zmiana korelacji między parą genów między datasetami |

### Dynamiczny
| Narzędzie | Opis |
|---|---|
| `execute_code` | Agent pisze własny JS gdy żadne narzędzie nie wystarczy |
