# Datierungsteam: Shiny-App für das Seminar „Datierungsmethoden“

Interaktive Lernumgebung zu ¹⁴C-Kalibration, Dendrochronologie, Münzen (TPQ), Schriftquellen (TAQ),
Stratigrafie und Bayes'scher Modellierung. Entstanden für ein Seminar an der Universität Tübingen.
Alle Beispieldaten sind **fiktiv**.

**App im Browser öffnen:** https://mrcbrnnr.github.io/datierung-seminar/

Die Web-Version läuft komplett im Browser (Shinylive, R als WebAssembly). Es muss weder R installiert
noch ein Server betrieben werden. Beim ersten Öffnen dauert das Laden etwa 20 bis 40 Sekunden.
Empfohlen: Laptop mit aktuellem Chrome, Edge, Firefox oder Safari.

## Was die App kann

| Tab | Inhalt |
|---|---|
| 1 · Kalibration & Modell | Proben und Schichten eingeben oder ein Szenario laden; Modell mit Stratigrafie, TPQ, TAQ und Dendro-Fenstern; Annahmen einzeln ein- und ausschalten; Vergleich mit der Ablagerung nur aus den ¹⁴C-Daten |
| 2 · Wiggle-Matching | Simulierte Ringsequenz an die IntCal20-Kurve legen, χ²-Test, Passungskurve |
| 3 · SPD & Kombination | Summenkurve (verschiedene Ereignisse) oder Kombination (dasselbe Ereignis, gewichtetes Mittel mit χ²-Test) |
| 4 · OxCal-Code | Erzeugt OxCal-Code aus dem Modell zum Vergleich |
| 5 · Altholz-Check | Wie viel älter ist eine Probe als die jüngste Probe ihrer Schicht (oder als die Ablagerung)? |
| 6 · Messwert → Alter | F14C in ¹⁴C-Alter umrechnen (und zurück) und kalibrieren |
| 7 · Dendro-Crossdating | Probe gegen eine **simulierte** Referenzchronologie verschieben (r, t-Wert, Gleichläufigkeit), Fälljahr-Fenster aus Waldkante, Splint oder Kernholz |

Szenarien in Tab 1: „Kastell Musterberg“ (mit Konflikt in Schicht C und mit gelöstem Konflikt),
„Plateau-Beispiel (Hallstatt-Zeit)“, „Altholz-Übung“ und „Kollektivgrab“ (Übung nach M. Hinz 2012).

Jahreszahlen dürfen in den Tabellen als Text eingegeben werden („480 v. Chr.“, „101 n. Chr.“).
Intern wird astronomisch gezählt wie in OxCal (Jahr 0 = 1 v. Chr.). Oben rechts lässt sich die
Anzeige zwischen v./n. Chr. und astronomischer Zählung umschalten.

## Repository

```
app/app.R                          die komplette App
.github/workflows/deploy.yml       baut die Web-Version und veröffentlicht sie über GitHub Pages
README.md                          diese Datei
```

Bei jeder Änderung an `app/app.R` auf dem Zweig `main` baut GitHub Actions die Seite automatisch neu
(einige Minuten). Voraussetzung: **Settings → Pages → Source: GitHub Actions**.
Die Kalibrationskurve IntCal20 wird beim Bau von intcal.org geladen und neben die App gelegt.

## Lokal starten

```r
install.packages("shiny")
shiny::runApp("app")
```

Die Kurve wird aus `app/intcal20.14c` gelesen. Fehlt die Datei, versucht die App lokal das Paket
`rcarbon` oder einen Download von intcal.org.

## Grenzen

- Das Modell ist ein didaktisches Monte-Carlo-Modell (Ziehen und Verwerfen), kein Ersatz für OxCal oder
  andere MCMC-Software. Ergebnisse sind ähnlich, aber nicht identisch.
- Nicht enthalten: Reservoireffekte, Marine20 und SHCal20, Ausreißermodelle (nur im OxCal-Export).
- Die Referenzchronologie in Tab 7 ist simuliert und ersetzt keine echte Standardchronologie.
- Die ¹⁴C-Alter der Fallstudie sind illustrativ gewählt.

## Quellen und Dank

- IntCal20: Reimer et al. 2020, *Radiocarbon* 62, 725–757. Bitte die Nutzungshinweise auf intcal.org beachten.
- Splintholz-Spannen (Eiche) in Tab 7 nach Tegel et al. 2022, *Frontiers in Ecology and Evolution* 10.
- Bayes'sche Grundlagen und OxCal-Bausteine: Bronk Ramsey 2009, *Radiocarbon* 51.
- Ideen und Übungsdaten (Kollektivgrab): Unterlagen von Martin Hinz zum Seminar „Absolute Chronologie und
  Isotopenforschung“ (2012), verwendet mit Quellenangabe.
