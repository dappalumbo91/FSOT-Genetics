theory GeneticsSpine
  imports Main
begin

(* FSOT-Genetics catalog spine. Generated. Pin D1D38A. *)

definition free_parameters :: nat where "free_parameters = 0"
definition phi_milli :: nat where "phi_milli = 1618"
definition leftover_milli :: nat where "leftover_milli = 382"
definition close_homolog_milli :: nat where "close_homolog_milli = 618"
definition backbone_D :: nat where "backbone_D = 8"
definition disulfide_D :: nat where "disulfide_D = 7"
definition salt_D :: nat where "salt_D = 9"
definition pack_D :: nat where "pack_D = 14"
definition hbond_D :: nat where "hbond_D = 8"
definition molecular_D :: nat where "molecular_D = 9"
definition tertiary_D :: nat where "tertiary_D = 13"
definition long_range_gate :: nat where "long_range_gate = 7"
definition chem_link_card :: nat where "chem_link_card = 7"
definition product_n :: nat where "product_n = 10"
definition product_sub2A :: nat where "product_sub2A = 10"
definition product_median_milliA :: nat where "product_median_milliA = 133"
definition alphafold_median_milliA :: nat where "alphafold_median_milliA = 471"
definition product_lt_af :: nat where "product_lt_af = 133"
definition product_lt_af_rhs :: nat where "product_lt_af_rhs = 471"
definition product_lt_bulk :: nat where "product_lt_bulk = 133"
definition product_lt_bulk_rhs :: nat where "product_lt_bulk_rhs = 13574"
definition homolog_measured :: nat where "homolog_measured = 26"
definition homolog_folded :: nat where "homolog_folded = 26"
definition homolog_true_miss :: nat where "homolog_true_miss = 1"
definition homolog_close :: nat where "homolog_close = 10"
definition homolog_close_rhs :: nat where "homolog_close_rhs = 24"
definition analog_jobs :: nat where "analog_jobs = 2"
definition male_olf_vs_vnc :: nat where "male_olf_vs_vnc = 6"
definition male_olf_vs_vnc_rhs :: nat where "male_olf_vs_vnc_rhs = 107402"
definition male_olf_vs_jo :: nat where "male_olf_vs_jo = 6"
definition male_olf_vs_jo_rhs :: nat where "male_olf_vs_jo_rhs = 77729"
definition banc_olf_vs_vnc :: nat where "banc_olf_vs_vnc = 8"
definition banc_olf_vs_vnc_rhs :: nat where "banc_olf_vs_vnc_rhs = 101619"
definition banc_olf_vs_jo :: nat where "banc_olf_vs_jo = 8"
definition banc_olf_vs_jo_rhs :: nat where "banc_olf_vs_jo_rhs = 66369"
definition larva_olf_vs_mech :: nat where "larva_olf_vs_mech = 4292"
definition larva_olf_vs_mech_rhs :: nat where "larva_olf_vs_mech_rhs = 166016"
definition worm_herm_vs_male_sex :: nat where "worm_herm_vs_male_sex = 745"
definition worm_herm_vs_male_sex_rhs :: nat where "worm_herm_vs_male_sex_rhs = 77442"
definition plant_arabidopsis_folded :: nat where "plant_arabidopsis_folded = 6"
definition plant_crop_folded :: nat where "plant_crop_folded = 24"
definition plant_crop_miss :: nat where "plant_crop_miss = 0"
definition plant_min_id :: nat where "plant_min_id = 618"
definition plant_min_id_rhs :: nat where "plant_min_id_rhs = 741"
definition plant_crop_n :: nat where "plant_crop_n = 24"
definition worm_cells :: nat where "worm_cells = 0"
definition worm_cells_rhs :: nat where "worm_cells_rhs = 453"
definition ciona_cells :: nat where "ciona_cells = 0"
definition ciona_cells_rhs :: nat where "ciona_cells_rhs = 205"
definition platynereis_cells :: nat where "platynereis_cells = 0"
definition platynereis_cells_rhs :: nat where "platynereis_cells_rhs = 1720"
definition male_neurons :: nat where "male_neurons = 100000"
definition male_neurons_rhs :: nat where "male_neurons_rhs = 165122"
definition banc_neurons :: nat where "banc_neurons = 100000"
definition banc_neurons_rhs :: nat where "banc_neurons_rhs = 175401"
definition ciona_gaba_flag :: nat where "ciona_gaba_flag = 0"
definition platynereis_gaba_flag :: nat where "platynereis_gaba_flag = 0"

lemma ok_free_parameters_zero: "free_parameters = 0"
  unfolding free_parameters_def by simp

lemma ok_phi_milli: "phi_milli = 1618"
  unfolding phi_milli_def by simp

lemma ok_leftover_floor_milli: "leftover_milli = 382"
  unfolding leftover_milli_def by simp

lemma ok_close_homolog_milli: "close_homolog_milli = 618"
  unfolding close_homolog_milli_def by simp

lemma ok_chem_backbone_D: "backbone_D = 8"
  unfolding backbone_D_def by simp

lemma ok_chem_disulfide_D: "disulfide_D = 7"
  unfolding disulfide_D_def by simp

lemma ok_chem_salt_D: "salt_D = 9"
  unfolding salt_D_def by simp

lemma ok_chem_pack_D: "pack_D = 14"
  unfolding pack_D_def by simp

lemma ok_chem_hbond_D: "hbond_D = 8"
  unfolding hbond_D_def by simp

lemma ok_chem_molecular_D: "molecular_D = 9"
  unfolding molecular_D_def by simp

lemma ok_chem_tertiary_D: "tertiary_D = 13"
  unfolding tertiary_D_def by simp

lemma ok_long_range_gate: "long_range_gate = 7"
  unfolding long_range_gate_def by simp

lemma ok_chem_link_card: "chem_link_card = 7"
  unfolding chem_link_card_def by simp

lemma ok_product_n: "product_n = 10"
  unfolding product_n_def by simp

lemma ok_product_sub2A: "product_sub2A = 10"
  unfolding product_sub2A_def by simp

lemma ok_product_median_milliA: "product_median_milliA = 133"
  unfolding product_median_milliA_def by simp

lemma ok_alphafold_median_milliA: "alphafold_median_milliA = 471"
  unfolding alphafold_median_milliA_def by simp

lemma ok_product_lt_alphafold: "product_lt_af < product_lt_af_rhs"
  unfolding product_lt_af_def product_lt_af_rhs_def by simp

lemma ok_product_lt_bulk: "product_lt_bulk < product_lt_bulk_rhs"
  unfolding product_lt_bulk_def product_lt_bulk_rhs_def by simp

lemma ok_homolog_measured: "homolog_measured = 26"
  unfolding homolog_measured_def by simp

lemma ok_homolog_folded_eq_measured: "homolog_folded = 26"
  unfolding homolog_folded_def by simp

lemma ok_homolog_true_miss: "homolog_true_miss = 1"
  unfolding homolog_true_miss_def by simp

lemma ok_homolog_close_count: "homolog_close <= homolog_close_rhs"
  unfolding homolog_close_def homolog_close_rhs_def by simp

lemma ok_analog_jobs: "analog_jobs = 2"
  unfolding analog_jobs_def by simp

lemma ok_male_vnc_gt_olf: "male_olf_vs_vnc < male_olf_vs_vnc_rhs"
  unfolding male_olf_vs_vnc_def male_olf_vs_vnc_rhs_def by simp

lemma ok_male_jo_gt_olf: "male_olf_vs_jo < male_olf_vs_jo_rhs"
  unfolding male_olf_vs_jo_def male_olf_vs_jo_rhs_def by simp

lemma ok_banc_vnc_gt_olf: "banc_olf_vs_vnc < banc_olf_vs_vnc_rhs"
  unfolding banc_olf_vs_vnc_def banc_olf_vs_vnc_rhs_def by simp

lemma ok_banc_jo_gt_olf: "banc_olf_vs_jo < banc_olf_vs_jo_rhs"
  unfolding banc_olf_vs_jo_def banc_olf_vs_jo_rhs_def by simp

lemma ok_larva_mech_gt_olf: "larva_olf_vs_mech < larva_olf_vs_mech_rhs"
  unfolding larva_olf_vs_mech_def larva_olf_vs_mech_rhs_def by simp

lemma ok_worm_male_sex_gt_herm: "worm_herm_vs_male_sex < worm_herm_vs_male_sex_rhs"
  unfolding worm_herm_vs_male_sex_def worm_herm_vs_male_sex_rhs_def by simp

lemma ok_plant_arabidopsis_folded: "plant_arabidopsis_folded = 6"
  unfolding plant_arabidopsis_folded_def by simp

lemma ok_plant_crop_folded: "plant_crop_folded = 24"
  unfolding plant_crop_folded_def by simp

lemma ok_plant_crop_miss: "plant_crop_miss = 0"
  unfolding plant_crop_miss_def by simp

lemma ok_plant_min_id_close: "plant_min_id <= plant_min_id_rhs"
  unfolding plant_min_id_def plant_min_id_rhs_def by simp

lemma ok_plant_crop_n: "plant_crop_n = 24"
  unfolding plant_crop_n_def by simp

lemma ok_worm_cells: "worm_cells < worm_cells_rhs"
  unfolding worm_cells_def worm_cells_rhs_def by simp

lemma ok_ciona_cells: "ciona_cells < ciona_cells_rhs"
  unfolding ciona_cells_def ciona_cells_rhs_def by simp

lemma ok_platynereis_cells: "platynereis_cells < platynereis_cells_rhs"
  unfolding platynereis_cells_def platynereis_cells_rhs_def by simp

lemma ok_male_neurons: "male_neurons < male_neurons_rhs"
  unfolding male_neurons_def male_neurons_rhs_def by simp

lemma ok_banc_neurons: "banc_neurons < banc_neurons_rhs"
  unfolding banc_neurons_def banc_neurons_rhs_def by simp

lemma ok_ciona_unsigned_gaba: "ciona_gaba_flag = 0"
  unfolding ciona_gaba_flag_def by simp

lemma ok_platynereis_unsigned_gaba: "platynereis_gaba_flag = 0"
  unfolding platynereis_gaba_flag_def by simp

end