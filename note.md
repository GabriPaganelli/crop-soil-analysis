# Note esplorative e considerazioni modellistiche

## Struttura generale del dataset

- Misure ripetute lungo il profilo verticale:
  
  \[
  \text{osservazioni nested in Field}
  \]

- Possibile struttura reale:

  \[
  \text{Field} \subset \text{Farm}
  \]

  ma `Farm` non è osservato.

- Ogni `Field` appartiene a una sola `Landuse`.

- Le osservazioni entro `Field` mostrano dipendenza intra-gruppo:
  - ICC moderato;
  - random intercept importante;
  - possibile random slope in profondità.

- La struttura verticale sembra essere il driver dominante del dataset:
  - forte decadimento con profondità;
  - eterogeneità tra field ancora elevata anche dopo aver corretto per depth.

---

# Profondità (`Bottom`)

## Questioni aperte

- `Bottom` rappresenta:
  - il limite inferiore dello strato?
  - oppure il valore medio dell’intervallo?

- Gli strati non sembrano avere ampiezza uniforme:
  - 20, 30, 40, 50, 60, 80 cm;
  - possibile perdita del 70/80.

- Se gli intervalli hanno spessori diversi:
  - la profondità potrebbe NON essere trattabile come asse continuo uniforme.

---

## Pattern osservati

- Forte decadimento verticale di SOC.

- Il decadimento NON sembra perfettamente lineare:
  - forte caduta iniziale;
  - possibile plateau in profondità.

- Il modello con `factor(Bottom)` suggerisce:
  - concavità del profilo;
  - dinamica più intensa nei primi 40–50 cm.

- Possibile eteroschedasticità lungo la profondità:
  
  \[
  Var(Y \mid Bottom)
  \downarrow
  \]

  all’aumentare della profondità.

---

## Random slopes

Dal mixed model:

\[
\text{slope} \propto \text{intercept}
\]

con correlazione negativa intercept–slope:

- i field con SOC superficiale elevato decadono più rapidamente;
- i field poveri mostrano profili più piatti.

Osservazione biologicamente plausibile e potenzialmente molto importante.

---

# Forme funzionali candidate

Da confrontare:

## Lineare

\[
Y \sim Bottom
\]

oppure

\[
\log(Y) \sim Bottom
\]

---

## Potenza

\[
\log(Y) \sim \log(Bottom)
\]

---

## Decadimento esponenziale

\[
Y \sim e^{-k \cdot Bottom}
\]

oppure con parametri dipendenti da covariate:

\[
Y \sim A(x)e^{-k(x)\cdot Bottom}
\]

---

# Struttura di correlazione

Possibili approcci:

- indipendenza condizionata dato `Field`;
- random intercept;
- random slope;
- correlazione seriale verticale:
  - AR(1),
  - exponential correlation,
  - ecc.

Possibile utilizzo di `nlme` per:
- strutture di correlazione;
- modellazione della varianza (`varIdent`, `varExp`, ...).

---

# Landuse e management

## Struttura osservata

`Landuse` è quasi deterministicamente funzione di:

\[
Landuse = f(Irrigate, Fertilised, OnFarm)
\]

tranne:
- `Landuse 3` e `Landuse 4`,
  che condividono la stessa tripla:

\[
(Irrigate,Fertilised,OnFarm)=(0,0,0)
\]

---

## Conseguenze

- `Landuse` e le tre variabili binarie NON sono indipendenti.
- Rischio:
  - collinearità;
  - singolarità;
  - non identificabilità di alcune interazioni.

In particolare:

\[
OnFarm = 1 \Rightarrow Fertilised = 1
\]

quindi:
- alcuni effetti non sono separabili.

---

## Aspetti interessanti

`Landuse 3` e `Landuse 4`:
- stessa gestione;
- texture diversa.

Possibile quasi-esperimento naturale per distinguere:
- effetto management;
- effetto texture.

---

# Texture

- Texture composizionale:
  
  \[
  Sand + Silt + Clay = 100
  \]

- Necessario evitare pseudoreplicazione composizionale:
  - ILR;
  - CLR.

- Forte correlazione tra texture e landuse:
  - quasi lineare nel triangle plot;
  - overlap parziale ma non completo.

Nel simplesso:
- `Landuse 3,6,7` concentrati verso un estremo;
- `5,2,1,4` distribuiti più centralmente.

Quindi:
- texture e management sono fortemente associati;
- ma non completamente confusi.

---

# Correlazioni tra variabili risposta

Correlazioni globali:
- SOC e nitrogeno: elevate;
- fosforo molto meno correlato.

Ma:

\[
Corr(Y_1,Y_2 \mid group) \neq Corr(Y_1,Y_2)
\]

Le correlazioni cambiano molto dopo stratificazione:
- per landuse;
- per field;
- per profondità.

Possibile:
- Simpson’s paradox;
- struttura gerarchica forte.

---

# Interpretazione causale

Dataset osservazionale.

Non è possibile distinguere chiaramente:

> “c’è più SOC perché è coltivato”

da

> “è coltivato perché il suolo ha già più SOC”.

Possibili:
- confondimento;
- selezione dei siti;
- bias ambientale.

Interpretazioni causali da formulare con cautela.

---

# Approcci modellistici candidati

1. LM semplici

2. Linear Mixed Models

3. Modelli con struttura di correlazione verticale

4. Depth non lineare:
   - spline;
   - polinomi;
   - funzioni esponenziali;
   - GAM/GAMM.

5. Decomposizione di `Landuse`:
   
   \[
   Irrigate * Fertilised * OnFarm
   \]

   con attenzione a:
   - singolarità;
   - celle vuote;
   - design sbilanciato.

6. Texture composizionale:
   - ILR;
   - CLR.

7. Modelli multivariati:
   - SOC,
   - nitrogeno,
   - fosforo congiuntamente.

8. Bayesian hierarchical models:
   - pooling parziale;
   - random slopes;
   - propagazione dell’incertezza;
   - strutture gerarchiche complesse.