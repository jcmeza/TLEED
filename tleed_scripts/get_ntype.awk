awk '
	/optimization parameters / {getline
	print
}' $* 
	
