using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random
using ITensorMPS
using Revise

ψ = cmps.CAMPS(5)
CampsPP.apply!(ψ, 2, P"XIZYY", π/2)
cmps.expectation(ψ, P"IIYII")
ITensorMPS.expect(ψ.mps,"Y")
ψ.Cdag