all:
	g++ -O3 -std=c++11 -o compiled/eclat_lamp_cmh_2k_bool -c baselines/eclat_lamp_cmh_2k.cpp
	g++ -O3 -std=c++11 -o compiled/eclat_lamp_cmh_Mk_bool -c baselines/eclat_lamp_cmh_Mk.cpp
	g++ -O3 -std=c++11 -o compiled/eclat_bonf_cmh_bool -c baselines/eclat_bonf_cmh.cpp
	g++ -O3 -std=c++11 -o compiled/eclat_lamp_cmh_np_bool -c baselines/eclat_lamp_cmh_np.cpp
	g++ -O3 -std=c++11 -o compiled/eclat_lamp_fisher_bool -c lamp/eclat_lamp_fisher.cpp
	g++ -O3 -std=c++11 -o compiled/eclat_lamp_chi_bool -c lamp/eclat_lamp_chi.cpp
	g++ -O3 -std=c++11 -o compiled/eclat_lamp_cmh_bool cmh/eclat_lamp_cmh.cpp
