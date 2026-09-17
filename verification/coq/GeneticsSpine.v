(* FSOT-Genetics catalog spine. Generated. Pin D1D38A. Prelude only. *)

Definition free_parameters := 0.
Definition phi_milli := 1618.
Definition leftover_milli := 382.
Definition close_homolog_milli := 618.
Definition backbone_D := 8.
Definition disulfide_D := 7.
Definition salt_D := 9.
Definition pack_D := 14.
Definition hbond_D := 8.
Definition molecular_D := 9.
Definition tertiary_D := 13.
Definition long_range_gate := 7.
Definition chem_link_card := 7.
Definition product_n := 10.
Definition product_sub2A := 10.
Definition product_median_milliA := 133.
Definition alphafold_median_milliA := 471.
Definition product_lt_af := 133.
Definition product_lt_af_rhs := 471.
Definition product_lt_bulk := 133.
Definition product_lt_bulk_rhs := 13574.
Definition homolog_measured := 26.
Definition homolog_folded := 26.
Definition homolog_true_miss := 1.
Definition homolog_close := 10.
Definition homolog_close_rhs := 24.
Definition analog_jobs := 2.
Definition male_olf_vs_vnc := 6.
Definition male_olf_vs_vnc_rhs := 107402.
Definition male_olf_vs_jo := 6.
Definition male_olf_vs_jo_rhs := 77729.
Definition banc_olf_vs_vnc := 8.
Definition banc_olf_vs_vnc_rhs := 101619.
Definition banc_olf_vs_jo := 8.
Definition banc_olf_vs_jo_rhs := 66369.
Definition larva_olf_vs_mech := 4292.
Definition larva_olf_vs_mech_rhs := 166016.
Definition worm_herm_vs_male_sex := 745.
Definition worm_herm_vs_male_sex_rhs := 77442.
Definition plant_arabidopsis_folded := 6.
Definition plant_crop_folded := 24.
Definition plant_crop_miss := 0.
Definition plant_min_id := 618.
Definition plant_min_id_rhs := 741.
Definition plant_crop_n := 24.
Definition worm_cells := 0.
Definition worm_cells_rhs := 453.
Definition ciona_cells := 0.
Definition ciona_cells_rhs := 205.
Definition platynereis_cells := 0.
Definition platynereis_cells_rhs := 1720.
Definition male_neurons := 100000.
Definition male_neurons_rhs := 165122.
Definition banc_neurons := 100000.
Definition banc_neurons_rhs := 175401.
Definition ciona_gaba_flag := 0.
Definition platynereis_gaba_flag := 0.

Lemma ok_free_parameters_zero : free_parameters = 0.
Proof. reflexivity. Qed.

Lemma ok_phi_milli : phi_milli = 1618.
Proof. reflexivity. Qed.

Lemma ok_leftover_floor_milli : leftover_milli = 382.
Proof. reflexivity. Qed.

Lemma ok_close_homolog_milli : close_homolog_milli = 618.
Proof. reflexivity. Qed.

Lemma ok_chem_backbone_D : backbone_D = 8.
Proof. reflexivity. Qed.

Lemma ok_chem_disulfide_D : disulfide_D = 7.
Proof. reflexivity. Qed.

Lemma ok_chem_salt_D : salt_D = 9.
Proof. reflexivity. Qed.

Lemma ok_chem_pack_D : pack_D = 14.
Proof. reflexivity. Qed.

Lemma ok_chem_hbond_D : hbond_D = 8.
Proof. reflexivity. Qed.

Lemma ok_chem_molecular_D : molecular_D = 9.
Proof. reflexivity. Qed.

Lemma ok_chem_tertiary_D : tertiary_D = 13.
Proof. reflexivity. Qed.

Lemma ok_long_range_gate : long_range_gate = 7.
Proof. reflexivity. Qed.

Lemma ok_chem_link_card : chem_link_card = 7.
Proof. reflexivity. Qed.

Lemma ok_product_n : product_n = 10.
Proof. reflexivity. Qed.

Lemma ok_product_sub2A : product_sub2A = 10.
Proof. reflexivity. Qed.

Lemma ok_product_median_milliA : product_median_milliA = 133.
Proof. reflexivity. Qed.

Lemma ok_alphafold_median_milliA : alphafold_median_milliA = 471.
Proof. reflexivity. Qed.

Lemma ok_product_lt_alphafold : Nat.ltb product_lt_af product_lt_af_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_product_lt_bulk : Nat.ltb product_lt_bulk product_lt_bulk_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_homolog_measured : homolog_measured = 26.
Proof. reflexivity. Qed.

Lemma ok_homolog_folded_eq_measured : homolog_folded = 26.
Proof. reflexivity. Qed.

Lemma ok_homolog_true_miss : homolog_true_miss = 1.
Proof. reflexivity. Qed.

Lemma ok_homolog_close_count : Nat.leb homolog_close homolog_close_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_analog_jobs : analog_jobs = 2.
Proof. reflexivity. Qed.

Lemma ok_male_vnc_gt_olf : Nat.ltb male_olf_vs_vnc male_olf_vs_vnc_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_male_jo_gt_olf : Nat.ltb male_olf_vs_jo male_olf_vs_jo_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_banc_vnc_gt_olf : Nat.ltb banc_olf_vs_vnc banc_olf_vs_vnc_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_banc_jo_gt_olf : Nat.ltb banc_olf_vs_jo banc_olf_vs_jo_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_larva_mech_gt_olf : Nat.ltb larva_olf_vs_mech larva_olf_vs_mech_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_worm_male_sex_gt_herm : Nat.ltb worm_herm_vs_male_sex worm_herm_vs_male_sex_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_plant_arabidopsis_folded : plant_arabidopsis_folded = 6.
Proof. reflexivity. Qed.

Lemma ok_plant_crop_folded : plant_crop_folded = 24.
Proof. reflexivity. Qed.

Lemma ok_plant_crop_miss : plant_crop_miss = 0.
Proof. reflexivity. Qed.

Lemma ok_plant_min_id_close : Nat.leb plant_min_id plant_min_id_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_plant_crop_n : plant_crop_n = 24.
Proof. reflexivity. Qed.

Lemma ok_worm_cells : Nat.ltb worm_cells worm_cells_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_ciona_cells : Nat.ltb ciona_cells ciona_cells_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_platynereis_cells : Nat.ltb platynereis_cells platynereis_cells_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_male_neurons : Nat.ltb male_neurons male_neurons_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_banc_neurons : Nat.ltb banc_neurons banc_neurons_rhs = true.
Proof. reflexivity. Qed.

Lemma ok_ciona_unsigned_gaba : ciona_gaba_flag = 0.
Proof. reflexivity. Qed.

Lemma ok_platynereis_unsigned_gaba : platynereis_gaba_flag = 0.
Proof. reflexivity. Qed.
