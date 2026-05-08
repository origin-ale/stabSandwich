clear
clear

# args: P χs t
# julia --threads=auto --project=. Sandwich_debug.jl XZ 1,2 2
# julia --project=. Sandwich_debug.jl ZXX 1,2,4
# echo -------------------------------------------------------------------
# julia --project=. Sandwich_debug.jl ZXXX 1,2,4,8
julia --project=. Sandwich_debug.jl ZZZZZZZ 1,2,4,8,16
# julia --project=. Sandwich_debug_deftime.jl ZZZZZZZ 14

# julia --threads=auto --project=. Sandwich_tomography.jl

# julia --threads=auto --project=. Debug_disentangler.jl