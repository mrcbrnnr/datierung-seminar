# Datierungsteam: Shiny-App für das Seminar „Datierungsmethoden“

Didaktische App zu ¹⁴C-Kalibration, Dendro-Fenstern, Münzen (TPQ), Schriftquellen (TAQ),
Stratigrafie, Altholz, Wiggle-Matching, Summenkurven und OxCal-Export.
Die Beispieldaten sind fiktiv.

**Live-Version:** `https://mrcbrnnr.github.io/datierung-seminar/`  
(Link nach der Einrichtung unten eintragen.)

Die Web-Version läuft komplett im Browser (Shinylive/WebAssembly). Es muss weder R installiert
noch ein Server betrieben werden. Beim ersten Öffnen dauert das Laden etwa eine halbe Minute.

## Einrichtung (einmalig)

1. Neues **öffentliches** Repository `datierung-seminar` auf GitHub anlegen und diese Dateien hochladen
   (`app/app.R`, `.github/workflows/deploy.yml`, `README.md`, `.gitignore`).
2. Im Repository: **Settings → Pages → Source: GitHub Actions** wählen.
3. Unter **Actions** den Lauf „Shiny-App als Webseite veröffentlichen“ abwarten (einige Minuten).
   Beim ersten Mal ggf. **Run workflow** drücken.
4. Die Adresse steht unter **Settings → Pages** und im Lauf unter „deploy“.

Bei jeder Änderung an `app/app.R` wird die Seite automatisch neu gebaut.

## Lokal starten (ohne Web-Version)

```r
install.packages("shiny")
shiny::runApp("app")
```

Die Kalibrationskurve wird aus `app/intcal20.14c` gelesen. Fehlt die Datei, versucht die App
lokal `rcarbon` oder einen Download von intcal.org. Die Web-Version bekommt die Datei im
Workflow automatisch.

## Quellen und Hinweise

- IntCal20: Reimer et al. 2020, Radiocarbon 62, 725–757. Bitte die Nutzungshinweise auf intcal.org beachten.
- Das Modell ist ein didaktisches Monte-Carlo-Modell und kein Ersatz für OxCal.
