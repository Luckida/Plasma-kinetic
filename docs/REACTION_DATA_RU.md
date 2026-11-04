# Зарегистрированные химические данные

Это все декларативные записи банка, включая альтернативы и оценки. Рабочая сеть конкретного запуска зависит от состава, флагов, доступных состояний и политики коэффициентов. Это не список всех возможных реакций природы.

Всего зарегистрировано **60 записей** в 5 группах. Источники и Q/области — в полных записях и CSV; таблица не назначает отсутствующие коэффициенты.

## H

Записей: 10. Выбор по умолчанию: `{}`. Участники-селекторы раскрыты в `selectors` этого JSON.

| ID и строка JSON | Реагенты → продукты | Тип | Коэффициенты | Статус / набор |
|---|---|---|---|---|
| [`UMIST22:74`](../data/reactions/H/reactions.json#L6) | H- + H → H2 + e | `heavy_chemistry` | [k(T), строка 2](../data/reactions/H/rate_coefficients.csv#L2) | запись банка |
| [`UMIST22:500`](../data/reactions/H/reactions.json#L22) | H + H2+ → H2 + H+ | `charge_exchange` | [k(T), строка 4](../data/reactions/H/rate_coefficients.csv#L4) | запись банка |
| [`UMIST22:1524`](../data/reactions/H/reactions.json#L38) | H2+ + e → 2 H | `dissociative_recombination` | [k(T), строка 5](../data/reactions/H/rate_coefficients.csv#L5) | запись банка |
| [`UMIST22:1577`](../data/reactions/H/reactions.json#L53) | H3+ + e → H2 + H | `dissociative_recombination` | [k(T), строка 6](../data/reactions/H/rate_coefficients.csv#L6) | запись банка |
| [`UMIST22:1578`](../data/reactions/H/reactions.json#L69) | H3+ + e → 3 H | `dissociative_recombination` | [k(T), строка 7](../data/reactions/H/rate_coefficients.csv#L7) | запись банка |
| [`UMIST22:3196`](../data/reactions/H/reactions.json#L84) | H2+ + H2 → H3+ + H | `heavy_chemistry` | [k(T), строка 8](../data/reactions/H/rate_coefficients.csv#L8) | запись банка |
| [`UMIST22:6353`](../data/reactions/H/reactions.json#L100) | H- + H+ → 2 H | `mutual_neutralization` | [k(T), строка 9](../data/reactions/H/rate_coefficients.csv#L9) | запись банка |
| [`UMIST22:6354`](../data/reactions/H/reactions.json#L115) | H- + H2+ → H2 + H | `mutual_neutralization` | [k(T), строка 10](../data/reactions/H/rate_coefficients.csv#L10) | запись банка |
| [`UMIST22:6357`](../data/reactions/H/reactions.json#L131) | H- + H3+ → 2 H + H2 | `mutual_neutralization` | [k(T), строка 11](../data/reactions/H/rate_coefficients.csv#L11) | запись банка |
| [`LS05:Hplus_H:CX`](../data/reactions/H/reactions.json#L147) | H+ + H → H + H+ | `charge_exchange` | [σ(E), строка 2](../data/reactions/H/cross_sections.csv#L2) | запись банка |

## H--He

Записей: 18. Выбор по умолчанию: `{"He_H_chemiionization": "RD82", "HeH_recombination": "UMIST22", "HeH_H_exchange": "UMIST22"}`. Участники-селекторы раскрыты в `selectors` этого JSON.

| ID и строка JSON | Реагенты → продукты | Тип | Коэффициенты | Статус / набор |
|---|---|---|---|---|
| [`UMIST22:469`](../data/reactions/H--He/reactions.json#L6) | H2 + He+ → He + H2+ | `charge_exchange` | [k(T), строка 2](../data/reactions/H--He/rate_coefficients.csv#L2) | запись банка |
| [`UMIST22:1711`](../data/reactions/H--He/reactions.json#L22) | HeH+ + e → He + H | `dissociative_recombination` | [k(T), строка 3](../data/reactions/H--He/rate_coefficients.csv#L3) | запись банка; UMIST22 |
| [`UMIST22:3272`](../data/reactions/H--He/reactions.json#L40) | H2 + He+ → H+ + H + He | `heavy_chemistry` | [k(T), строка 4](../data/reactions/H--He/rate_coefficients.csv#L4) | запись банка |
| [`UMIST22:3273`](../data/reactions/H--He/reactions.json#L57) | H2 + HeH+ → He + H3+ | `heavy_chemistry` | [k(T), строка 5](../data/reactions/H--He/rate_coefficients.csv#L5) | запись банка |
| [`UMIST22:3682`](../data/reactions/H--He/reactions.json#L73) | H + HeH+ → He + H2+ | `heavy_chemistry` | [k(T), строка 6](../data/reactions/H--He/rate_coefficients.csv#L6) | запись банка; UMIST22 |
| [`UMIST22:6364`](../data/reactions/H--He/reactions.json#L91) | H- + He+ → H + He | `mutual_neutralization` | [k(T), строка 7](../data/reactions/H--He/rate_coefficients.csv#L7) | запись банка |
| [`RD82:He3S_H:PI`](../data/reactions/H--He/reactions.json#L107) | He(2^3S) + H → He + H+ + e | `penning` | [k(T), строка 8](../data/reactions/H--He/rate_coefficients.csv#L8) | запись банка; RD82 |
| [`RD82:He3S_H:AI`](../data/reactions/H--He/reactions.json#L126) | He(2^3S) + H → HeH+ + e | `associative_ionization` | [k(T), строка 9](../data/reactions/H--He/rate_coefficients.csv#L9) | запись банка; RD82 |
| [`GM25:He3S_H:PI`](../data/reactions/H--He/reactions.json#L144) | He(2^3S) + H → He + H+ + e | `penning` | [k(T), строка 10](../data/reactions/H--He/rate_coefficients.csv#L10) | запись банка; GM25 |
| [`GM25:He3S_H:AI`](../data/reactions/H--He/reactions.json#L163) | He(2^3S) + H → HeH+ + e | `associative_ionization` | [k(T), строка 14](../data/reactions/H--He/rate_coefficients.csv#L14) | запись банка; GM25 |
| [`GM25:He1S_H:PI`](../data/reactions/H--He/reactions.json#L181) | He(2^1S) + H → He + H+ + e | `penning` | [k(T), строка 18](../data/reactions/H--He/rate_coefficients.csv#L18) | запись банка; GM25 |
| [`GM25:He1S_H:AI`](../data/reactions/H--He/reactions.json#L200) | He(2^1S) + H → HeH+ + e | `associative_ionization` | [k(T), строка 22](../data/reactions/H--He/rate_coefficients.csv#L22) | запись банка; GM25 |
| [`GM25:He3S_Hplus:CT`](../data/reactions/H--He/reactions.json#L218) | He(2^3S) + H+ → He+ + H | `charge_exchange` | [k(T), строка 26](../data/reactions/H--He/rate_coefficients.csv#L26) | оценка, по умолчанию выключена |
| [`GM25:He1S_Hplus:CT`](../data/reactions/H--He/reactions.json#L235) | He(2^1S) + H+ → He+ + H | `charge_exchange` | [k(T), строка 30](../data/reactions/H--He/rate_coefficients.csv#L30) | оценка, по умолчанию выключена |
| [`GM25:Heplus_H:RCT`](../data/reactions/H--He/reactions.json#L252) | He+ + H → He + H+ | `charge_exchange` | [k(T), строка 34](../data/reactions/H--He/rate_coefficients.csv#L34) | запись банка; с фотоном |
| [`GM25:HeH_e:DR`](../data/reactions/H--He/reactions.json#L269) | HeH+ + e → He + H | `dissociative_recombination` | [k(T), строка 38](../data/reactions/H--He/rate_coefficients.csv#L38) | запись банка; GM25 |
| [`GM25:HeH_H:exchange`](../data/reactions/H--He/reactions.json#L287) | HeH+ + H → He + H2+ | `heavy_chemistry` | [k(T), строка 42](../data/reactions/H--He/rate_coefficients.csv#L42) | запись банка; GM25 |
| [`GM25:H2plus_He:exchange`](../data/reactions/H--He/reactions.json#L305) | H2+ + He → HeH+ + H | `heavy_chemistry` | [k(T), строка 46](../data/reactions/H--He/rate_coefficients.csv#L46) | запись банка |

## H--Li

Записей: 3. Выбор по умолчанию: `{}`. Участники-селекторы раскрыты в `selectors` этого JSON.

| ID и строка JSON | Реагенты → продукты | Тип | Коэффициенты | Статус / набор |
|---|---|---|---|---|
| [`GP98:Li_Hplus:CT`](../data/reactions/H--Li/reactions.json#L10) | Li + H+ → Li+ + H | `charge_exchange` | [k(T), строка 2](../data/reactions/H--Li/rate_coefficients.csv#L2) | запись банка |
| [`GP98:Li_Hplus:RCT`](../data/reactions/H--Li/reactions.json#L27) | Li + H+ → Li+ + H | `charge_exchange` | [k(T), строка 3](../data/reactions/H--Li/rate_coefficients.csv#L3) | запись банка; с фотоном |
| [`GP98:H_Liplus:inverse_CT`](../data/reactions/H--Li/reactions.json#L44) | Li+ + H → Li + H+ | `charge_exchange` | [k(T), строка 2](../data/reactions/H--Li/rate_coefficients.csv#L2) | запись банка |

## He

Записей: 4. Выбор по умолчанию: `{}`. Участники-селекторы раскрыты в `selectors` этого JSON.

| ID и строка JSON | Реагенты → продукты | Тип | Коэффициенты | Статус / набор |
|---|---|---|---|---|
| [`HE25:15:PI`](../data/reactions/He/reactions.json#L18) | 2 He(2^3S) → He + He+ + e | `penning` | [k(T), строка 2](../data/reactions/He/rate_coefficients.csv#L2) | запись банка |
| [`HE25:14:AI`](../data/reactions/He/reactions.json#L34) | 2 He(2^3S) → He2+ + e | `associative_ionization` | [k(T), строка 3](../data/reactions/He/rate_coefficients.csv#L3) | запись банка |
| [`HE25:9:DR_ground`](../data/reactions/He/reactions.json#L49) | He2+ + e → 2 He | `dissociative_recombination` | [k(T), строка 4](../data/reactions/He/rate_coefficients.csv#L4) | запись банка |
| [`HE25:10:DR_triplet`](../data/reactions/He/reactions.json#L64) | He2+ + e → He(2^3S) + He | `dissociative_recombination` | [k(T), строка 5](../data/reactions/He/rate_coefficients.csv#L5) | запись банка |

## He--Li

Записей: 25. Выбор по умолчанию: `{}`. Участники-селекторы раскрыты в `selectors` этого JSON.

| ID и строка JSON | Реагенты → продукты | Тип | Коэффициенты | Статус / набор |
|---|---|---|---|---|
| [`BR15:2->3`](../data/reactions/He--Li/reactions.json#L54) | Li(2s) + He+ → Li+ + He(2^3S) | `charge_exchange` | [σ(E), строка 2](../data/reactions/He--Li/cross_sections.csv#L2); [k(T), строка 2](../data/reactions/He--Li/rate_coefficients.csv#L2) | запись банка |
| [`BR15:3->2`](../data/reactions/He--Li/reactions.json#L71) | Li+ + He(2^3S) → Li(2s) + He+ | `charge_exchange` | [σ(E), строка 312](../data/reactions/He--Li/cross_sections.csv#L312); [k(T), строка 407](../data/reactions/He--Li/rate_coefficients.csv#L407) | запись банка |
| [`BR15:2->4`](../data/reactions/He--Li/reactions.json#L88) | Li(2s) + He+ → Li+ + He(2^1S) | `charge_exchange` | [σ(E), строка 621](../data/reactions/He--Li/cross_sections.csv#L621); [k(T), строка 83](../data/reactions/He--Li/rate_coefficients.csv#L83) | запись банка |
| [`BR15:4->2`](../data/reactions/He--Li/reactions.json#L105) | Li+ + He(2^1S) → Li(2s) + He+ | `charge_exchange` | [σ(E), строка 1594](../data/reactions/He--Li/cross_sections.csv#L1594); [k(T), строка 650](../data/reactions/He--Li/rate_coefficients.csv#L650) | запись банка |
| [`BR15:2->5`](../data/reactions/He--Li/reactions.json#L122) | Li(2s) + He+ → Li+ + He(2^3P) | `charge_exchange` | [σ(E), строка 2566](../data/reactions/He--Li/cross_sections.csv#L2566); [k(T), строка 164](../data/reactions/He--Li/rate_coefficients.csv#L164) | запись банка |
| [`BR15:5->2`](../data/reactions/He--Li/reactions.json#L139) | Li+ + He(2^3P) → Li(2s) + He+ | `charge_exchange` | [σ(E), строка 2892](../data/reactions/He--Li/cross_sections.csv#L2892); [k(T), строка 893](../data/reactions/He--Li/rate_coefficients.csv#L893) | запись банка |
| [`BR15:2->6`](../data/reactions/He--Li/reactions.json#L156) | Li(2s) + He+ → Li(2p) + He+ | `heavy_quenching` | [σ(E), строка 3217](../data/reactions/He--Li/cross_sections.csv#L3217); [k(T), строка 245](../data/reactions/He--Li/rate_coefficients.csv#L245) | запись банка |
| [`BR15:6->2`](../data/reactions/He--Li/reactions.json#L173) | Li(2p) + He+ → Li(2s) + He+ | `heavy_quenching` | [σ(E), строка 4463](../data/reactions/He--Li/cross_sections.csv#L4463); [k(T), строка 1136](../data/reactions/He--Li/rate_coefficients.csv#L1136) | запись банка |
| [`BR15:2->7`](../data/reactions/He--Li/reactions.json#L190) | Li(2s) + He+ → Li+ + He(2^1P) | `charge_exchange` | [σ(E), строка 5708](../data/reactions/He--Li/cross_sections.csv#L5708); [k(T), строка 326](../data/reactions/He--Li/rate_coefficients.csv#L326) | запись банка |
| [`BR15:7->2`](../data/reactions/He--Li/reactions.json#L207) | Li+ + He(2^1P) → Li(2s) + He+ | `charge_exchange` | [σ(E), строка 6056](../data/reactions/He--Li/cross_sections.csv#L6056); [k(T), строка 1541](../data/reactions/He--Li/rate_coefficients.csv#L1541) | запись банка |
| [`BR15:3->5`](../data/reactions/He--Li/reactions.json#L224) | Li+ + He(2^3S) → Li+ + He(2^3P) | `heavy_quenching` | [σ(E), строка 6403](../data/reactions/He--Li/cross_sections.csv#L6403); [k(T), строка 488](../data/reactions/He--Li/rate_coefficients.csv#L488) | запись банка |
| [`BR15:5->3`](../data/reactions/He--Li/reactions.json#L241) | Li+ + He(2^3P) → Li+ + He(2^3S) | `heavy_quenching` | [σ(E), строка 6642](../data/reactions/He--Li/cross_sections.csv#L6642); [k(T), строка 974](../data/reactions/He--Li/rate_coefficients.csv#L974) | запись банка |
| [`BR15:3->6`](../data/reactions/He--Li/reactions.json#L258) | Li+ + He(2^3S) → Li(2p) + He+ | `charge_exchange` | [σ(E), строка 6880](../data/reactions/He--Li/cross_sections.csv#L6880); [k(T), строка 569](../data/reactions/He--Li/rate_coefficients.csv#L569) | запись банка |
| [`BR15:6->3`](../data/reactions/He--Li/reactions.json#L275) | Li(2p) + He+ → Li+ + He(2^3S) | `charge_exchange` | [σ(E), строка 7161](../data/reactions/He--Li/cross_sections.csv#L7161); [k(T), строка 1217](../data/reactions/He--Li/rate_coefficients.csv#L1217) | запись банка |
| [`BR15:4->6`](../data/reactions/He--Li/reactions.json#L292) | Li+ + He(2^1S) → Li(2p) + He+ | `charge_exchange` | [σ(E), строка 7441](../data/reactions/He--Li/cross_sections.csv#L7441); [k(T), строка 731](../data/reactions/He--Li/rate_coefficients.csv#L731) | запись банка |
| [`BR15:6->4`](../data/reactions/He--Li/reactions.json#L309) | Li(2p) + He+ → Li+ + He(2^1S) | `charge_exchange` | [σ(E), строка 7940](../data/reactions/He--Li/cross_sections.csv#L7940); [k(T), строка 1298](../data/reactions/He--Li/rate_coefficients.csv#L1298) | запись банка |
| [`BR15:4->7`](../data/reactions/He--Li/reactions.json#L326) | Li+ + He(2^1S) → Li+ + He(2^1P) | `heavy_quenching` | [σ(E), строка 8438](../data/reactions/He--Li/cross_sections.csv#L8438); [k(T), строка 812](../data/reactions/He--Li/rate_coefficients.csv#L812) | запись банка |
| [`BR15:7->4`](../data/reactions/He--Li/reactions.json#L343) | Li+ + He(2^1P) → Li+ + He(2^1S) | `heavy_quenching` | [σ(E), строка 9660](../data/reactions/He--Li/cross_sections.csv#L9660); [k(T), строка 1622](../data/reactions/He--Li/rate_coefficients.csv#L1622) | запись банка |
| [`BR15:5->6`](../data/reactions/He--Li/reactions.json#L360) | Li+ + He(2^3P) → Li(2p) + He+ | `charge_exchange` | [σ(E), строка 10881](../data/reactions/He--Li/cross_sections.csv#L10881); [k(T), строка 1055](../data/reactions/He--Li/rate_coefficients.csv#L1055) | запись банка |
| [`BR15:6->5`](../data/reactions/He--Li/reactions.json#L377) | Li(2p) + He+ → Li+ + He(2^3P) | `charge_exchange` | [σ(E), строка 11225](../data/reactions/He--Li/cross_sections.csv#L11225); [k(T), строка 1379](../data/reactions/He--Li/rate_coefficients.csv#L1379) | запись банка |
| [`BR15:6->7`](../data/reactions/He--Li/reactions.json#L394) | Li(2p) + He+ → Li+ + He(2^1P) | `charge_exchange` | [σ(E), строка 11568](../data/reactions/He--Li/cross_sections.csv#L11568); [k(T), строка 1460](../data/reactions/He--Li/rate_coefficients.csv#L1460) | запись банка |
| [`BR15:7->6`](../data/reactions/He--Li/reactions.json#L411) | Li+ + He(2^1P) → Li(2p) + He+ | `charge_exchange` | [σ(E), строка 12788](../data/reactions/He--Li/cross_sections.csv#L12788); [k(T), строка 1703](../data/reactions/He--Li/rate_coefficients.csv#L1703) | запись банка |
| [`S23:He3S_Li:PI_model`](../data/reactions/He--Li/reactions.json#L428) | He(2^3S) + Li(2s) → He + Li+ + e | `penning` | [σ(E), строка 14007](../data/reactions/He--Li/cross_sections.csv#L14007) | оценка, по умолчанию выключена |
| [`S23:He1S_Li:PI_model`](../data/reactions/He--Li/reactions.json#L446) | He(2^1S) + Li(2s) → He + Li+ + e | `penning` | [σ(E), строка 14248](../data/reactions/He--Li/cross_sections.csv#L14248) | оценка, по умолчанию выключена |
| [`S23:He3P_Li:PI_model`](../data/reactions/He--Li/reactions.json#L464) | He(2^3P) + Li(2s) → He + Li+ + e | `penning` | [σ(E), строка 14489](../data/reactions/He--Li/cross_sections.csv#L14489) | оценка, по умолчанию выключена |
