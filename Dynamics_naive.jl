using Revise

using CampsPP
import PauliPropagation as pp
import CliffordMPS as cmps

N = 101
χ = 128

ψ = cmps.CAMPS(N)
obs = pp.PauliSum(pp.PauliString(N, [:Z], [N÷2]))

output = "output/NaiveDynamics.txt"
open(output, "w") do f
    println(f, "# N=$N χ=$χ obs=Z[N÷2]")
end

ψ, k, s = camps_rndrotation_dynamics(ψ, χ, obs, output; showprogress=true)