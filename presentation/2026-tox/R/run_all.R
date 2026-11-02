# Rebuild every data file and figure for the 2026-tox talk.
#   cd presentation/2026-tox && Rscript R/run_all.R
# then assemble the deck:  python3 build_deck.py
for (f in c("s1_library.R", "s2_ms.R", "s3_ici.R", "s4_obesity.R",
            "s5_doac.R", "s6_merigolix.R", "s7_diagrams.R")) {
  message("== ", f)
  status <- system2("Rscript", file.path("R", f))
  if (status != 0) stop(f, " failed")
}
