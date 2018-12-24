do for[a=1:246]{
	pl 'rosen_func.txt' u 1:2:3 every :::a::a w p pt 7 ps 5 palette
	pause 0.1
}
pause -1
