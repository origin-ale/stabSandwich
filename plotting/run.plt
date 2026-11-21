# Shared run selection for the gnuplot scripts (logic in run_dirs.sh).
#
# Set `datadir` (the folder the simulation script writes, e.g.
# "output/Comptime_prob_tm") before loading. Optional override: run=<run id>
# (gnuplot -e "run='20261120-160717-Comptime_prob_tm'" …); an explicit `dir`
# (or dir_campspp/dir_pp) takes precedence over it.
#
# Defines
#   rundir       default input directory: output/runs/<run>/files/<datadir>, the
#                fetched run's own copy, for the run given as `run`, else for the
#                newest fetched run that wrote datadir; datadir itself if no
#                fetched run did
#   plot_dir(d)  output/plots/<run id>, the plot folder for the input dir(s) d
#                (space-separated), created if missing. For dirs outside
#                output/runs the "run id" is the dir's own name, e.g.
#                output/old/carlos_varseed -> output/plots/carlos_varseed.
if (exists("run")) {
  rundir = system(sprintf("plotting/run_dirs.sh data '%s' '%s'", datadir, run))
} else {
  rundir = system(sprintf("plotting/run_dirs.sh data '%s'", datadir))
}
#   method_pos(f, m)  method m's 0-based position in the methods= header of
#                file f, -1 if the run didn't include it, 0 for a file without
#                the header (old data, one method per file)
if (strlen(rundir) == 0) { exit status 1 }
plot_dir(d) = system("plotting/run_dirs.sh plotdir ".d)
method_pos(f, m) = system(sprintf("plotting/run_dirs.sh method '%s' '%s'", f, m)) + 0
