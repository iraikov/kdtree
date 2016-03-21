kdtree: kdtree.mlb kdtree.sml
	mlton -cc-opt "-O3" -mlb-path-var "TENSOR_LIB $(HOME)/src/ML/tensor" -mlb-path-var "ARRAYSORT_LIB $(HOME)/src/ML/arraysort"  $< 
kdtreetest: kdtreetest.mlb kdtree.sml kdtreetest.sml
	mlton -const 'Exn.keepHistory true' -cc-opt "-O3" -mlb-path-var "TENSOR_LIB $(HOME)/src/ML/tensor" -mlb-path-var "ARRAYSORT_LIB $(HOME)/src/ML/arraysort"  $< 



