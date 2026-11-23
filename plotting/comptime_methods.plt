# Shared method selection for plot_comptime_prob_tm.plt and
# plot_comptime_prob_tm_box.plt. Load after run.plt.
#
# Methods: campspp, camps, pp, in this order in the key and, in the box plots,
# from left to right at each n. A method the run didn't include is left out.
#
# dir_<m>  directory holding method m's comptime_prob_tm.txt (and
#          gctime_prob_tm.txt). Default: rundir for every method. Once any
#          dir_<m> is given, only the methods given one are plotted, so old
#          per-method folders don't get mixed with the newest run.
# j_<m>    method m's 0-based position in its file's `methods`; -1 leaves it
#          out (default: from the file's methods= header, -1 if the run didn't
#          include m, 0 for a file without the header).
#
# Defines mlist (the methods to plot, space-separated), dirof(m), jof(m),
# pick(m, a, b, c) (a for campspp, b for camps, c for pp) and plotdir; exits
# with an error if no method is left.
all_methods = "campspp camps pp"
dirs_given = 0
do for [m in all_methods] { if (exists("dir_".m)) { dirs_given = 1 } }
do for [m in all_methods] {
  if (!exists("dir_".m)) { eval sprintf("dir_%s = '%s'", m, dirs_given ? "" : rundir) }
  if (!exists("j_".m)) {
    if (value("dir_".m) eq "") {
      eval sprintf("j_%s = -1", m)
    } else {
      eval sprintf("j_%s = %d", m, method_pos(value("dir_".m).'/comptime_prob_tm.txt', m))
    }
  }
}
mlist = ""
plotdirs = ""
do for [m in all_methods] {
  if (value("j_".m) >= 0) { mlist = mlist.(mlist eq "" ? "" : " ").m; plotdirs = plotdirs." ".value("dir_".m) }
}
if (words(mlist) == 0) {
  print "comptime_methods.plt: none of campspp, camps, pp in the run's methods"
  exit status 1
}
dirof(m) = value("dir_".m)
jof(m) = value("j_".m)
pick(m, a, b, c) = (m eq "campspp") ? a : ((m eq "camps") ? b : c)
plotdir = plot_dir(plotdirs)
