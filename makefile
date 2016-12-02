all:
	g++ -O3 -std=c++11 -o compiled/facs_2k baselines/eclat_lamp_cmh_2k.cpp
	g++ -O3 -std=c++11 -o compiled/facs_Mk baselines/eclat_lamp_cmh_Mk.cpp
	g++ -O3 -std=c++11 -o compiled/bonf_cmh baselines/eclat_bonf_cmh.cpp
	g++ -O3 -std=c++11 -o compiled/lamp_chi lamp/eclat_lamp_chi.cpp
	g++ -O3 -std=c++11 -o compiled/facs cmh/eclat_lamp_cmh.cpp
