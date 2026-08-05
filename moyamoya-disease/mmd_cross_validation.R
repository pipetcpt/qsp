## mmd_cross_validation.R -- step 2 of the two-implementation check.
## Run `python3 mmd_cross_validation.py` first: it writes /tmp/py_ref.json.
## This script re-computes the same quantities in mrgsolve and prints the
## paired comparison.  Every equation in the model is implemented twice,
## independently, and this is the file that holds them to each other.

suppressMessages({library(mrgsolve); library(dplyr); library(jsonlite)})
src <- readLines("mmd_mrgsolve_model.R")
k <- grep("^if \\(identical\\(environment", src)
eval(parse(text = paste(src[1:(k-1)], collapse="\n")))
py <- fromJSON("/tmp/py_ref.json")
o <- run_mmd("adult_isch", days=3650, delta=5)
h <- run_mmd("adult_haem", days=3650, delta=5)
b <- run_mmd("adult_isch", days=3285+90, delta=0.05,
             surgery=list(kind="direct", day=3285), pars=list(G_LEAK=0.35))
MAPK <- c(CBFA="CBFA",CBFF="CBFF",CBFWS_="CBFWS",OEFA="OEFA",CVR_INTR="CVRINT",
          CVR_MEAS="CVRMEAS",INFA="INFA",ISCH="ISCH",STEN="STEN",gS="gS",
          GMOYA="GMOYA",GPVA="GPVA",ANEU="ANEU",HEMH="HEMH",REMOD="REMOD",
          SIG_PVA="SIGPVA")
runs <- list(adult_isch=o, adult_haem=h)
res <- data.frame()
for (nm in names(py)) {
  parts <- strsplit(nm, "\\|")[[1]]
  if (parts[1] == "bypass") next
  d <- as.numeric(parts[2]); key <- parts[3]
  df <- runs[[parts[1]]]
  got <- approx(df$time, df[[MAPK[[key]]]], d)$y
  want <- py[[nm]]
  rel <- if (abs(want) > 1e-6) abs(got-want)/abs(want)*100 else abs(got-want)
  res <- rbind(res, data.frame(item=nm, R=got, py=want, pct=rel))
}
m <- b$time >= 3285
bp <- c(peak_HYPER_REL=max(b$HYPERREL[m]), peak_CBFF=max(b$CBFF[m]),
        d90_CBFWS=approx(b$time,b$CBFWS,3285+90)$y,
        d90_GBYP=approx(b$time,b$GBYP,3285+90)$y)
for (nn in names(bp)) {
  want <- py[[paste0("bypass|",nn)]]
  res <- rbind(res, data.frame(item=paste0("bypass|",nn), R=bp[[nn]], py=want,
                               pct=abs(bp[[nn]]-want)/max(abs(want),1e-9)*100))
}
res$pct <- round(res$pct,2)
res$R <- round(res$R,4); res$py <- round(res$py,4)
cat("\n== mrgsolve vs Python reference:", nrow(res), "paired values ==\n")
print(res[order(-res$pct),][1:18,], row.names=FALSE)
cat(sprintf("\n  median |diff| %.3f %%   mean %.3f %%   max %.3f %%\n",
            median(res$pct), mean(res$pct), max(res$pct)))
cat(sprintf("  within 1%%: %d/%d   within 5%%: %d/%d\n",
            sum(res$pct<1), nrow(res), sum(res$pct<5), nrow(res)))
