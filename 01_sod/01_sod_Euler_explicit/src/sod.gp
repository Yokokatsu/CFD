# interactive.gp
set xlabel "x [m]"
set ylabel "value"

set xrange [0:1]
set yrange [0:1]
set key right top
set grid

nsteps = 300
dt = 0.001

do for [i=1:nsteps] {
    t = i * dt
    set title sprintf("shock tube %04d", i)
    filename = sprintf("./../result/sod_%04d.dat", i)
    plot filename using 1:2 with lines lw 2 linecolor rgb "red" title "velocity u [m/s]", \
         filename using 1:3 with lines lw 2 linecolor rgb "blue" title "density ρ [kg/m³]", \
         filename using 1:4 with lines lw 2 linecolor rgb "green" title "pressure p [Pa]"
    pause 0.01
}

pause -1