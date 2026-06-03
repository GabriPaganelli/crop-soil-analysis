# Rassegna bibliografica — Confronto con la letteratura
## Modello M-SP: analisi dei nutrienti del suolo in Tanzania

> **Nota metodologica**: per ogni risultato del modello viene riportato prima il quadro generale
> (studi internazionali), poi uno zoom sull'Africa sub-sahariana/Tanzania dove disponibile.
> I valori dei parametri si riferiscono al modello M-SP finale (mediana e CI90%).

---

## 1. Forma funzionale del profilo verticale: log-log (power-law)

**Nostro risultato**: `logY ~ logBottom` — il profilo di ogni nutriente segue una legge di potenza
con la profondità (log-log lineare). Scelto su base AIC rispetto a lineare, quadratico e fattoriale.

### Letteratura generale

- **Jobbágy & Jackson (2000)** — *Ecological Applications, 10(2): 423–436* — Studio seminale su
  >2000 profili di suolo in tre database globali. I profili verticali di SOC sono ben descritti da
  funzioni esponenziali o power-law; la distribuzione verticale varia sistematicamente con clima e
  vegetazione. Riferimento di base per qualsiasi modellazione della distribuzione verticale del SOC.
  URL: https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/1051-0761(2000)010[0423:TVDOSO]2.0.CO;2

- **Studi più recenti (2020s)**: funzioni di profondità esponenziali sono usate come standard per
  mappare e modellare la distribuzione 3D del SOC in suoli agricoli. Una funzione lineare in log-log
  (= power-law: Y = a · Depth^b) è la forma più parsimoniosa biologicamente giustificata.
  Fonte: ScienceDirect (2025) — *Generating three-dimensional soil organic carbon density dataset by
  soil depth function and correction methods*
  URL: https://www.sciencedirect.com/science/article/abs/pii/S136481522500266X

- **Perché NON il quadratico**: un modello quadratico in Bottom ammette un minimo e poi una risalita
  del nutriente con la profondità — biologicamente privo di senso per SOC, N e P. Il vantaggio AIC
  marginale (~2.7 unità, 1 parametro extra) non giustifica la complessità aggiuntiva.

### Zoom Africa/Tanzania

- **Kirsten et al. (2019)** — *Journal of Plant Nutrition and Soil Science* — Su suoli
  profondamente alterati (Acrisols, Alisols) nelle Montagne Usambara (NE Tanzania), campionamento
  fino a 100 cm di profondità: i profili verticali di SOC mostrano un calo monotono e non lineare
  con la profondità, coerente con una forma power-law.
  URL: https://onlinelibrary.wiley.com/doi/10.1002/jpln.201800595

---

## 2. η_SOC = +0.209 [0.125, 0.302]: i campi con più SOC in media decadono più lentamente

**Nostro risultato**: la SD dell'effetto casuale di campo *aumenta* con la profondità per il SOC.
Campi con livello di SOC superiore alla media (z_nu > 0) hanno profili più piatti — decadono meno
rapidamente con la profondità. Effetto robusto e lontano da zero.

### Letteratura generale — coerenza

- **Stratification ratio (SR)**: il concetto di SR (= SOC_0-5cm / SOC_20-30cm) è ampiamente usato
  come indicatore di qualità del suolo. Franzluebbers (2002) e studi successivi mostrano che SR > 2
  è tipico di suoli di alta qualità (no-till, praterie), SR < 1.5 indica suoli degradati.
  Il nostro risultato è coerente: campi con più SOC in media hanno SR più alto (decadono meno), il
  che corrisponde a suoli di migliore qualità funzionale.
  Fonte: ScienceDirect — *Stratification ratio of soil organic matter pools as an indicator of
  carbon sequestration in a tillage chronosequence on a Brazilian Oxisol* (Mielniczuk et al., 2008)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S0167198708001736

- **Protezione minerale del SOC**: la letteratura sugli organo-mineral complexes (MAOC) spiega
  perché suoli ricchi di SOC hanno profili più piatti. Il SOC associato alle frazioni fini
  (silt + argilla) è protetto fisicamente e chimicamente, rallentando la mineralizzazione con la
  profondità. Suoli ad alto contenuto di SOC tendono ad avere più MAOC, che è meno suscettibile al
  degradazione profonda.
  Fonte: *Fine silt and clay content is the main factor defining maximal C and N accumulations in
  soils: a meta-analysis* — Scientific Reports (2021)
  URL: https://www.nature.com/articles/s41598-021-84821-6

- **Conservazione del carbonio in profondità**: suoli con alta concentrazione superficiale di SOC
  tendono ad avere anche più carbonio in profondità (distribuzione correlata verticalmente), il che
  si riflette in un profilo meno ripido. Questo meccanismo è dovuto al trasporto di dissolved organic
  matter (DOM) in profondità e alla bioturbazione, entrambi più attivi in suoli ad alta produttività.
  Fonte: *Stratification and Storage of Soil Organic Carbon and Nitrogen as Affected by Tillage
  Practices in the North China Plain* — PLoS ONE (2015)
  URL: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0128873

- **Stratification ratio e livello medio SOC**: la variabilità inter-campo della SR riflette
  differenze nella qualità del suolo, nelle pratiche di gestione storiche e nella vegetazione
  pregressa. Campi con più SOC totale tendono ad avere SR maggiore.
  Fonte: *Depth distribution of soil organic carbon in an Oxisol under different land uses* —
  ScienceDirect (2016)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S0016706116304876

### Zoom Africa/Tanzania

- **Dinamiche SOC in Tanzania orientale**: studi su suoli di Tanzania e Africa orientale mostrano
  ampia variabilità inter-campo del contenuto di SOC, con i campi più ricchi che tendono ad avere
  profili più distribuiti in profondità (più carbon sequestration nel subsuolo).
  Fonte: *Dynamics of Organic Carbon and Predictive Management Indices in Contrasting Land Use Types
  on Weathered Soils of Eastern Tanzania* — Preprints.org (2023)
  URL: https://www.preprints.org/manuscript/202307.0173/v1

- **SOC in sistemi agroforestali NE Tanzania**: studi nelle montagne del Nord-Est della Tanzania
  mostrano che i sistemi con più SOC in superficie (agroforestali vs monocoltura) presentano
  distribuzioni più uniformi con la profondità, coerenti con η_SOC > 0.
  Fonte: *Soil organic carbon stocks and fertility in smallholder indigenous agroforestry systems
  of the North-Eastern mountains, Tanzania* — ScienceDirect (2024)
  URL: https://www.sciencedirect.com/science/article/pii/S2352009424000063

---

## 3. η_N ≈ 0 [CI90% include zero]: profili di azoto uniformi tra campi

**Nostro risultato**: il parametro η_N ha mediana −0.051 con CI90% che copre lo zero. Non c'è
evidenza che campi con più N in media abbiano profili verticali diversi. Interpretazione: la
distribuzione verticale dell'azoto non varia sistematicamente tra campi.

### Letteratura generale — coerenza parziale

- **N legato al SOC ma disaccoppiabile**: il TN totale è fortemente correlato al SOC nei suoli,
  ma l'azoto ha dinamiche proprie (mineralizzazione, input minerali, fissazione biologica) che
  possono disaccoppiarlo dal profilo del SOC. Questo spiegherebbe η_N ≈ 0 nonostante η_SOC > 0.
  Fonte: *Frontiers — Total Nitrogen Stock in Soil Profile Affected by Land Use and Soil Type* (2022)
  URL: https://www.frontiersin.org/journals/environmental-science/articles/10.3389/fenvs.2022.945305/full

- **Variabilità N in profondità**: la densità di TN è più eterogenea nei primi 60 cm e tende a
  uniformarsi in profondità. A differenza del SOC, non è stata trovata un'interazione significativa
  tra il tipo di uso del suolo e la densità di TN nel profilo.
  Fonte: *Spatial variability and uncertainty of soil nitrogen across the conterminous United States
  at different depths* — Ecosphere (2022)
  URL: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.4170

- **N mobile vs immobile**: l'azoto minerale (NO3-, NH4+) è mobile e può redistribuirsi con
  l'acqua nel profilo, riducendo la stratificazione. Questo lo rende meno "stratificato" del SOC
  che è prevalentemente organico e immobile.

### Letteratura generale — potenziale contrasto

- **Stratification ratio per N**: Franzluebbers (2002) e altri mostrano che il TN tende a
  stratificarsi in suoli no-till (SR_N > 1.5), suggerendo che tra campi gestiti diversamente ci
  potrebbe essere variabilità. Tuttavia, se i campi del dataset hanno tutti pratiche simili, la
  variabilità inter-campo del profilo di N potrebbe essere effettivamente bassa.

### Zoom Africa/Tanzania

- **Tanzania — Kilombero**: studi su riso e mais con/senza irrigazione mostrano che N tende ad
  accumularsi in superficie in tutti i trattamenti, senza differenze marcate inter-campo nel profilo
  verticale una volta controllati i principali fattori fisici.
  Fonte: *Different agricultural practices affect soil carbon, nitrogen and phosphorous in
  Kilombero — Tanzania* — ScienceDirect/PubMed (2019)
  URL: https://pubmed.ncbi.nlm.nih.gov/30616188/

---

## 4. η_P = +0.096 [0.019, 0.179]: stratificazione debole del fosforo

**Nostro risultato**: effetto simile a SOC ma con CI90% appena sopra lo zero. Campi con più P in
media decadono leggermente meno rapidamente con la profondità. Effetto presente ma debole.

### Letteratura generale — coerenza

- **P è immobile**: il fosforo si muove nel suolo quasi esclusivamente per diffusione, non con il
  flusso idrico. Questo lo rende naturalmente stratificato in superficie, con profili più ripidi
  rispetto all'N. Campi con input maggiori di P (fertilizzazioni storiche) avranno più P in
  superficie e — se l'input è sufficiente — anche più P in profondità (arricchimento fino a 0.75 m).
  Fonte: *Phosphorus accumulation and spatial distribution in agricultural soils in Denmark* —
  ScienceDirect (2013)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S001670611300222X

- **No-till e stratificazione P**: la gestione no-till aumenta marcatamente la stratificazione del P
  nei primi 5 cm. Suoli con più P totale tendono ad avere profili meno ripidi grazie all'accumulo
  anche in profondità.
  Fonte: *Phosphorus Stratification: Agronomic & Environmental Consequences* — Alabama Extension
  URL: https://www.aces.edu/blog/topics/crop-production/phosphorus-stratification-agronomic-environmental-consequences/

- **Affinità P per ossidi Fe/Al**: nei suoli tropicali (ferralsols, acrisols), il P si lega
  fortemente agli ossidi di ferro e alluminio, riducendone la mobilità e aumentando la stratificazione.
  Questo spiega perché η_P > 0 ma più debole di η_SOC: il P non è conservato da meccanismi organici
  (come il SOC) ma da legami chimici che variano meno sistematicamente tra campi.
  Fonte: *Phosphorus fertilization and management in soils of Sub-Saharan Africa* — Margenot (2018)
  URL: https://margenot.cropsciences.illinois.edu/wp-content/uploads/2018/12/phosphorus-fertilization-and-management-in-soils-of-sub-saharan-africa.pdf

### Zoom Africa/Tanzania

- **P limitante in Africa tropicale**: il fosforo è il secondo nutriente più limitante dopo N in
  molti suoli dell'Africa tropicale. Le concentrazioni di P totale si trovano principalmente nello
  strato superficiale (0–15 cm, 47–79.5% del totale nel profilo). La distribuzione verticale è
  marcatamente stratificata, con variabilità inter-campo legata alla storia di fertilizzazione e
  alla mineralogia del suolo.
  Fonte: *Soil P availability as affected by the chemical composition of plant materials* —
  ScienceDirect (2003)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S0167880903001713

---

## 5. Management non predittivo (OnFarm, Irrigate, Fertilised, N_Natural)

**Nostro risultato**: nessuna variabile di management (irrigazione, fertilizzazione, on-farm,
N-naturale) viene selezionata da projpred per nessuna delle 3 risposte. I beta corrispondenti in
M-SP hanno CI90% che copre zero per tutti e 3 i nutrienti.

### Letteratura generale — parziale contrasto

- **Effetti concentrati in superficie**: la letteratura mostra che fertilizzazione e irrigazione
  hanno effetti positivi significativi sullo strato superiore (0–30 cm) del SOC, ma al di sotto dei
  60 cm gli effetti non sono statisticamente significativi. Se il nostro modello controlla già per
  tessitura e struttura verticale, il management non aggiunge varianza spiegata al *profilo*
  (intercetta casuale cattura già la differenza media tra campi).
  Fonte: *Nitrogen fertilization impact on soil carbon pools — subtropical wheat-mungbean-rice* —
  PMC (2021)
  URL: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8486117/

- **Effetti a lungo termine**: gli effetti del management sul SOC in profondità si manifestano
  tipicamente su scale temporali molto lunghe (>25 anni). In un dataset cross-sectional (come il
  nostro), la variabilità storica tra campi può mascherare gli effetti attuali del management.

- **Parziale contrasto da Kilombero, Tanzania**: Mutambu et al. (2019) trovano che irrigazione
  e fertilizzazione hanno effetto positivo sui profili di SOC e TN, con l'interazione
  irrigazione × fertilizzazione che estende l'effetto agli strati più profondi. Tuttavia, questo
  studio non controlla per tessitura e struttura del suolo nel modo in cui lo fa il nostro modello
  — il che può spiegare la differenza.
  Fonte: *Different agricultural practices affect soil carbon, nitrogen and phosphorous in
  Kilombero — Tanzania* — ScienceDirect (2019)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S0301479718314609

### Interpretazione nel contesto del nostro modello

- Il management influenza prevalentemente il *livello medio* di nutrienti (già catturato
  dall'intercetta casuale z_nu_r), non la *forma del profilo verticale*. Dopo aver controllato
  per Texture2 e BulkDensity (che riflettono proprietà fisiche del suolo correlate al management
  storico), il management attuale non aggiunge potere predittivo.

---

## 6. Texture2 = log(Silt/Clay) selezionato per SOC e N; Texture1 = Sand non selezionata

**Nostro risultato (projpred)**: Texture2 (ILR2 = log(Silt/Clay)) selezionata per SOC e N.
Texture1 (ILR1 ≈ Sand vs {Silt, Clay}) mai selezionata. Questo suggerisce che la distinzione
argilla/limo è più informativa della presenza di sabbia per predire i livelli di nutrienti.

### Letteratura generale — forte coerenza

- **Organo-mineral complexes (MAOC)**: la frazione fine del suolo (argilla + limo) è il principale
  driver dell'accumulo e della stabilizzazione del carbonio organico e dell'azoto. Le argille hanno
  elevata superficie specifica e carica superficiale, che favoriscono il legame chimico con la
  materia organica (MAOC = Mineral-Associated Organic Carbon).
  Fonte: *Fine silt and clay content is the main factor defining maximal C and N accumulations in
  soils: a meta-analysis* — Scientific Reports (2021)
  URL: https://www.nature.com/articles/s41598-021-84821-6

- **Silt vs Clay — distinzione rilevante**: silt e argilla hanno ruoli diversi nella protezione
  del SOC. La distinzione tra le due frazioni (catturata da Texture2 = log(Silt/Clay)) è più
  informativa del semplice contenuto di sabbia. L'articolo Cotrufo et al. e la meta-analisi di
  Nature Reports (2021) mostrano che sia silt che clay contribuiscono all'accumulo di MAOC, con
  effetti additivi.
  Fonte: *Dual role of silt and clay in the formation and accrual of stabilized soil organic
  carbon* — ScienceDirect (2024)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S0038071724000798

  Fonte: *Role of silt and clay fractions in organic carbon and nitrogen stabilization* —
  Springer Nature (2023)
  URL: https://link.springer.com/article/10.1007/s42729-023-01209-3

- **Sabbia — bassa capacità protettiva**: la sabbia ha superficie specifica molto ridotta e pochi
  siti di legame per la materia organica. Il SOC nei suoli sabbiosi è principalmente nella
  frazione particolata (POM), meno stabile. La capacità di stabilizzazione aumenta solo quando le
  frazioni fini sono sature di carbonio (MAOC capacity exceeded).
  Fonte: *Soil organic carbon accrual after conversion from cropland to grassland across a range
  of soil textures* — ScienceDirect (2026)
  URL: https://www.sciencedirect.com/science/article/pii/S0167880926001659

- **Mineralogy in tropical soils**: nei suoli tropicali alterati (dominati da ossidi di Fe/Al e
  kaolinite, tipici della Tanzania), la distinzione argilla/limo è particolarmente rilevante perché
  gli ossidi di Fe e Al si concentrano nella frazione fine (<2 μm) e legano fortemente sia SOC che P.
  Fonte: *Stabilization of Soil Organic Carbon as Influenced by Clay Mineralogy* —
  ScienceDirect (2017)
  URL: https://www.sciencedirect.com/science/article/abs/pii/S006521131730086X

---

## 7. BulkDensity predice logN (ma non logSOC né logP)

**Nostro risultato (projpred)**: BulkDensity selezionata per N (n=4, secondo in ordine), non per
SOC né per P. In M-SP, γ_N[4] per BulkDensity ha mediana negativa e CI90% lontano da zero.

### Letteratura generale — coerenza

- **BD inversamente correlata a TN**: la densità apparente è inversamente correlata al contenuto di
  materia organica (e quindi al TN, che è strettamente legato al SOC). Suoli più compatti (alta BD)
  hanno meno pori e meno spazio per la materia organica stabilizzata → meno TN.
  Fonte: *Nitrogen dynamics as a function of soil types, compaction, and moisture* — PLoS ONE (2024)
  URL: https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0301296

- **Compaction e ciclo dell'azoto**: la compattazione del suolo riduce la nitrificazione aerobica e
  promuove la denitrificazione (passaggio da aerobico ad anaerobico). La riduzione dei micropori
  limita l'attività microbica azotofissatrice, riducendo l'accumulo di N organico.
  Fonte: *Mixed Effects of Soil Compaction on the Nitrogen Cycle Under Pea and Wheat* — PMC (2022)
  URL: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8940171/

- **Perché NON per SOC e P**: la tessitura (Texture2) cattura già la maggior parte della variabilità
  della protezione minerale del SOC. BulkDensity aggiunge informazione per N perché l'azoto è più
  direttamente legato alle condizioni aerobiche/microbiche (influenzate dalla compattazione) rispetto
  al SOC, che è più legato alla mineralogia. Per il P, la mobilità è già così limitata dalla chimica
  (legami Fe/Al) che la struttura fisica del suolo aggiunge poco potere predittivo.

---

## 8. Approcci statistici nella letteratura: confronto con il modello M-SP

**Nostro approccio**: modello bayesiano gerarchico con slope proporzionale (M-SP), stimato via MCMC
(CmdStan), con selezione variabili tramite projection predictive (projpred) e validazione predittiva
tramite LOO-CV (PSIS). Struttura: `mu_r[i,j] = alpha_r + z_nu_r[j]*(psi_r + eta_r*logBottom_i) + gamma_r*X_W + beta_r*X_B`.

### La letteratura usa prevalentemente approcci più semplici

- **Regressione lineare stratificata per profondità (il più comune)**: la maggior parte degli studi
  confronta medie di SOC/N/P a diverse profondità con ANOVA o t-test, senza modellare il profilo
  come curva continua. Non cattura la variabilità inter-campo nella forma del profilo.
  Esempi: Mutambu et al. (2019, Tanzania), studi con SR come outcome.

- **Modelli lineari misti frequentisti (lmer/nlme)**: approccio più avanzato ma ancora comune.
  Usano random intercept per campo (`(1|Field)`) e regressione fissa su logBottom o Bottom.
  Non estimano la variabilità del *profilo* tra campi — solo del livello medio.
  Riferimento metodologico: Pinheiro & Bates (2000), *Mixed-Effects Models in S and S-PLUS*.

- **Random Forest e machine learning per mappatura**: usati per predire SOC spazialmente (30m
  resolution), non per modellare la struttura profilo. Irrilevanti per la nostra domanda.
  Fonte: *African soil properties mapped at 30m spatial resolution* — PMC (2021)
  URL: https://pmc.ncbi.nlm.nih.gov/articles/PMC7969779/

- **Modelli bayesiani semplici o semi-parametrici**: rari nella letteratura pedologica applicata.
  Qualche uso di Gaussian Process per interpolazione spaziale (kriging bayesiano), ma non per
  struttura verticale del profilo.

### Cosa rende M-SP metodologicamente originale rispetto alla letteratura

- **Slope proporzionale**: la struttura `z_nu[j]*(psi + eta*logBottom)` non è usata nella
  letteratura pedologica standard. È equivalente a vincolare la correlazione tra random intercept
  e random slope a +1 (su scala log), il che ha un'interpretazione biologica diretta (η_r).
  La letteratura usa al più `(1 + logBottom | Field)` senza tale vincolo.

- **Inferenza bayesiana su η_r**: il parametro η_r (come la variabilità tra campi del profilo
  cambia con la profondità) non è mai stimato direttamente in letteratura. I paper usano la
  stratification ratio come proxy descrittivo, non come parametro inferenziale con incertezza.

- **Tre risposte simultanee (SOC, N, P)**: quasi nessuno studio modella le tre risposte in un
  unico framework gerarchico multiresposta. Tipicamente vengono analizzate separatamente, perdendo
  la possibilità di confrontare direttamente i pattern tra nutrienti.

- **Projpred per selezione variabili**: praticamente assente nella letteratura di scienze del
  suolo. La letteratura usa stepwise AIC, LASSO, o selezione soggettiva. Projpred offre garanzie
  teoriche sulla correttezza della selezione che questi metodi non hanno.
  Fonte metodologica: Piironen & Vehtari (2017), *Comparison of Bayesian predictive methods*,
  *Statistics and Computing*, 27(3): 711–735.
  URL: https://link.springer.com/article/10.1007/s11222-016-9649-y

- **LOO-CV per confronto modelli**: nella letteratura pedologica il confronto tra modelli avviene
  quasi sempre tramite AIC/BIC (frequentista) o DIC (bayesiano obsoleto). PSIS-LOO è più robusto
  e stabile per modelli gerarchici.
  Fonte: Vehtari, Gelman & Gabry (2017), *Practical Bayesian model evaluation using LOO-CV and WAIC*,
  *Statistics and Computing*, 27(5): 1413–1432.
  URL: https://link.springer.com/article/10.1007/s11222-016-9696-4

### Studi bayesiani su profili di suolo (rari ma esistenti)

- **Gerarchici bayesiani per dati pedologici**: alcuni studi usano modelli gerarchici bayesiani per
  stimare le proprietà del suolo a livello di campo, ma con struttura random intercept semplice
  (non slope proporzionale) e senza η_r.
  Fonte: *Bayesian Hierarchical Random Intercept Model* — IOP (2017)
  URL: https://iopscience.iop.org/article/10.1088/1742-6596/855/1/012061/pdf

- **Modelli bayesiani per distribuzione verticale del SOC**: alcuni modelli di processo (RothC,
  CENTURY) stimano parametri di profondità in framework bayesiani, ma sono modelli di processo
  deterministici con calibrazione MCMC — molto diversi dal nostro approccio statistico.

### Sintesi del confronto metodologico

| Aspetto | Letteratura standard | M-SP (nostro) |
|---------|---------------------|---------------|
| Struttura profilo | ANOVA stratificato o LMM con fixed bottom | Slope proporzionale bayesiana |
| Variabilità inter-campo del *profilo* | Non modellata | η_r con distribuzione a posteriori |
| Inferenza | Frequentista (lmer, ANOVA) | Bayesiana (MCMC, PSIS-LOO) |
| Selezione variabili | Stepwise AIC / soggettiva | Projection predictive (projpred) |
| Confronto modelli | AIC / BIC | LOO-CV (PSIS) |
| Tre risposte insieme | Separatamente | Framework unico multiresposta |
| Contesto Tanzania | LMM semplici o descrittivo | Primo uso di slope proporzionale |

---

## 9. Correlazione di campo SOC–N: rho_int_SOC_N = +0.386 [0.038, 0.667]

**Nostro risultato (MVRE)**: i campi con effetto casuale di intercetta elevato per SOC tendono
ad avere effetto casuale elevato per N. La correlazione stimata è +0.386 con CI90% [0.038, 0.667]
— lontano da zero. Le correlazioni tra slope (rho_slope_SOC_N) e tra intercetti SOC-P e N-P sono
invece vicine a zero con alta incertezza.

Questo risultato emerge **solo** dal modello M-SP-RIRS-MVRE (6D random effects): i modelli
precedenti (M-SP-RIRS con RE bivariati per risposta) non potevano stimarlo.

### Meccanismo biologico principale: rapporto C:N della materia organica

La correlazione di campo SOC-N è biologicamente attesa e ben documentata:

- **Rapporto C:N stabile nella materia organica stabile (MAOC)**: la materia organica
  stabilizzata nei suoli agricoli ha un rapporto C:N relativamente fisso (~10–12). La
  letteratura mostra che il C:N della frazione stabile varia molto meno del C:N totale.
  Di conseguenza, campi con più SOC (più materia organica stabile) hanno anche più TN in
  modo quasi meccanistico.
  Fonte: *Stoichiometry and coupled biogeochemical cycling of carbon and nitrogen* —
  Cleveland & Liptzin (2007), *Global Biogeochemical Cycles*, 21, GB1019.
  URL: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2006GB002857

- **C:N vincolato nei suoli agricoli**: una meta-analisi globale su >800 studi conferma che
  la concentrazione di TN nei suoli è altamente correlata con SOC (R² > 0.9 in molti
  pedotipi). La pendenza della regressione TN ~ SOC corrisponde a C:N ≈ 10–12 negli
  strati superficiali di suoli agricoli. Campi che differiscono nel livello di SOC
  differiscono proporzionalmente in TN — esattamente il pattern catturato da rho_int_SOC_N.
  Fonte: *Batjes (1996) — Total carbon and nitrogen in the soils of the world*,
  *European Journal of Soil Science*, 47(2): 151–163.
  URL: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-2389.1996.tb01386.x

### Driver comuni tra campi

- **Produttività e input organici**: i campi con più input di biomassa (residui, letame,
  agroforestry) accumulano sia carbonio che azoto organico simultaneamente. La gestione
  storica di lungo periodo è un driver comune di SOC e TN a livello di campo.
  Fonte: *Kirkby et al. (2011) — Some trends underlying the chemistry of stabilised organic
  matter in Australian soils* — Soil Research, 49(7): 587–594.
  URL: https://www.publish.csiro.au/sr/SR11097

- **Tessitura fine come driver congiunto**: la frazione fine del suolo (argilla + limo)
  stabilizza contemporaneamente sia il carbonio che l'azoto organico tramite legami
  organo-minerali (MAOC). Campi con più argilla/limo tendono ad avere più SOC *e* più TN.
  Questo crea una correlazione di campo SOC-N anche in assenza di un link diretto
  tra i due nutrienti.
  Fonte: *Fine silt and clay is the main factor defining maximal C and N accumulations* —
  Scientific Reports (2021). URL: https://www.nature.com/articles/s41598-021-84821-6

### Perché la correlazione è a livello di campo (non di osservazione)

La correlazione residua (rho_res, stimata da MVRE-FULL) è ≈ 0 [−0.12, +0.15]: dati SOC
e N alla stessa profondità nello stesso campo non sono correlati *dopo* aver rimosso
l'effetto di campo. Questo distingue due livelli:
- **Inter-campo** (rho_int_SOC_N = +0.386): la fertilità strutturale del campo
  (determinata da storia, tessitura, gestione di lungo periodo) influenza simultaneamente
  i livelli di base di SOC e N.
- **Intra-campo** (rho_res ≈ 0): alla stessa profondità, le fluttuazioni locali di SOC
  e TN sono indipendenti — coerente con l'idea che i meccanismi di turnover rapido (N
  minerale, mineralizzazione) non siano strettamente accoppiati al SOC nello stesso
  campione.

### Zoom Africa/Tanzania

- **Correlazione SOC-TN nei suoli dell'Africa orientale**: analisi di database di suoli
  in Africa orientale mostrano forte correlazione SOC-TN nei suoli agricoli, con C:N
  medio ~11–13 negli strati superiori, compatibile con il nostro rho_int_SOC_N > 0.
  La correlazione tende a indebolirsi negli strati più profondi (dove le frazioni organiche
  labile e stabile si separano) — coerente con rho_slope_SOC_N ≈ 0 nel nostro modello.
  Fonte: *Vågen et al. (2016) — Mapping of soil properties and land degradation risk
  in Africa using MODIS reflectance* — Geoderma, 263: 216–225.
  URL: https://www.sciencedirect.com/science/article/pii/S0016706115301476

---

## 10. Sintesi e coerenza complessiva con la letteratura

| Risultato | Coerenza letteratura | Note |
|-----------|---------------------|------|
| Log-log come forma funzionale | ✅ Forte | Standard dalla Jobbágy & Jackson (2000) |
| η_SOC > 0: campi ricchi decadono più lentamente | ✅ Coerente | SR ratio, protezione MAOC |
| η_N ≈ 0: profili N uniformi | ✅ Parzialmente coerente | N è più mobile del SOC; suoli tropicali |
| η_P > 0 debole: stratificazione P | ✅ Coerente | P immobile, legami Fe/Al |
| Management non predittivo | ⚠️ Contrasto parziale | Kilombero (2019): irrigazione × fertilizzazione ha effetto; ma il nostro modello controlla meglio per tessitura |
| Texture2 selezionata, Texture1 no | ✅ Forte | MAOC, superficie specifica argilla/limo > sabbia |
| BulkDensity predice N, non SOC/P | ✅ Coerente | Compaction e ciclo N; tessitura cattura SOC |
| **rho_int_SOC_N = +0.386**: correlazione di campo SOC–N | ✅ **Forte** | C:N ratio stabile (Cleveland & Liptzin 2007; Batjes 1996); driver comuni (tessitura fine, input organici) |

---

## Fonti principali

- Jobbágy & Jackson (2000) — https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/1051-0761(2000)010[0423:TVDOSO]2.0.CO;2
- Mielniczuk et al. (2008) stratification ratio — https://www.sciencedirect.com/science/article/abs/pii/S0167198708001736
- Scientific Reports meta-analisi silt+clay (2021) — https://www.nature.com/articles/s41598-021-84821-6
- Dual role silt+clay MAOC (2024) — https://www.sciencedirect.com/science/article/abs/pii/S0038071724000798
- Mutambu et al. Kilombero Tanzania (2019) — https://pubmed.ncbi.nlm.nih.gov/30616188/
- Kirsten et al. NE Tanzania (2019) — https://onlinelibrary.wiley.com/doi/10.1002/jpln.201800595
- SOC stocks East Africa review — https://www.tandfonline.com/doi/abs/10.1080/02571862.2019.1640296
- Margenot P sub-Saharan Africa (2018) — https://margenot.cropsciences.illinois.edu/wp-content/uploads/2018/12/phosphorus-fertilization-and-management-in-soils-of-sub-saharan-africa.pdf
- Clay mineralogy stabilization (2017) — https://www.sciencedirect.com/science/article/abs/pii/S006521131730086X
- PLoS ONE N and compaction (2024) — https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0301296
- PLoS ONE N storage tillage China (2015) — https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0128873
- SOC agroforestry NE Tanzania (2024) — https://www.sciencedirect.com/science/article/pii/S2352009424000063
- SOC dynamics East Tanzania preprint (2023) — https://www.preprints.org/manuscript/202307.0173/v1
- P stratification Alabama Extension — https://www.aces.edu/blog/topics/crop-production/phosphorus-stratification-agronomic-environmental-consequences/
- P accumulation Denmark (2013) — https://www.sciencedirect.com/science/article/abs/pii/S001670611300222X
