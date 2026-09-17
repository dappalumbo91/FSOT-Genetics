module FSOTGenetics

(* Generated catalog. Pin D1D38A. Nat milliscale. *)

let free_parameters : nat = 0
let phi_milli : nat = 1618
let leftover_milli : nat = 382
let close_homolog_milli : nat = 618
let backbone_D : nat = 8
let disulfide_D : nat = 7
let salt_D : nat = 9
let pack_D : nat = 14
let hbond_D : nat = 8
let molecular_D : nat = 9
let tertiary_D : nat = 13
let long_range_gate : nat = 7
let chem_link_card : nat = 7
let product_n : nat = 10
let product_sub2A : nat = 10
let product_median_milliA : nat = 133
let alphafold_median_milliA : nat = 471
let product_lt_af : nat = 133
let product_lt_af_rhs : nat = 471
let product_lt_bulk : nat = 133
let product_lt_bulk_rhs : nat = 13574
let homolog_measured : nat = 26
let homolog_folded : nat = 26
let homolog_true_miss : nat = 1
let homolog_close : nat = 10
let homolog_close_rhs : nat = 24
let analog_jobs : nat = 2
let male_olf_vs_vnc : nat = 6
let male_olf_vs_vnc_rhs : nat = 107402
let male_olf_vs_jo : nat = 6
let male_olf_vs_jo_rhs : nat = 77729
let banc_olf_vs_vnc : nat = 8
let banc_olf_vs_vnc_rhs : nat = 101619
let banc_olf_vs_jo : nat = 8
let banc_olf_vs_jo_rhs : nat = 66369
let larva_olf_vs_mech : nat = 4292
let larva_olf_vs_mech_rhs : nat = 166016
let worm_herm_vs_male_sex : nat = 745
let worm_herm_vs_male_sex_rhs : nat = 77442
let plant_arabidopsis_folded : nat = 6
let plant_crop_folded : nat = 24
let plant_crop_miss : nat = 0
let plant_min_id : nat = 618
let plant_min_id_rhs : nat = 741
let plant_crop_n : nat = 24
let worm_cells : nat = 0
let worm_cells_rhs : nat = 453
let ciona_cells : nat = 0
let ciona_cells_rhs : nat = 205
let platynereis_cells : nat = 0
let platynereis_cells_rhs : nat = 1720
let male_neurons : nat = 100000
let male_neurons_rhs : nat = 165122
let banc_neurons : nat = 100000
let banc_neurons_rhs : nat = 175401
let ciona_gaba_flag : nat = 0
let platynereis_gaba_flag : nat = 0

let _ = assert (free_parameters = 0)
let _ = assert (phi_milli = 1618)
let _ = assert (leftover_milli = 382)
let _ = assert (close_homolog_milli = 618)
let _ = assert (backbone_D = 8)
let _ = assert (disulfide_D = 7)
let _ = assert (salt_D = 9)
let _ = assert (pack_D = 14)
let _ = assert (hbond_D = 8)
let _ = assert (molecular_D = 9)
let _ = assert (tertiary_D = 13)
let _ = assert (long_range_gate = 7)
let _ = assert (chem_link_card = 7)
let _ = assert (product_n = 10)
let _ = assert (product_sub2A = 10)
let _ = assert (product_median_milliA = 133)
let _ = assert (alphafold_median_milliA = 471)
let _ = assert (product_lt_af < product_lt_af_rhs)
let _ = assert (product_lt_bulk < product_lt_bulk_rhs)
let _ = assert (homolog_measured = 26)
let _ = assert (homolog_folded = 26)
let _ = assert (homolog_true_miss = 1)
let _ = assert (homolog_close <= homolog_close_rhs)
let _ = assert (analog_jobs = 2)
let _ = assert (male_olf_vs_vnc < male_olf_vs_vnc_rhs)
let _ = assert (male_olf_vs_jo < male_olf_vs_jo_rhs)
let _ = assert (banc_olf_vs_vnc < banc_olf_vs_vnc_rhs)
let _ = assert (banc_olf_vs_jo < banc_olf_vs_jo_rhs)
let _ = assert (larva_olf_vs_mech < larva_olf_vs_mech_rhs)
let _ = assert (worm_herm_vs_male_sex < worm_herm_vs_male_sex_rhs)
let _ = assert (plant_arabidopsis_folded = 6)
let _ = assert (plant_crop_folded = 24)
let _ = assert (plant_crop_miss = 0)
let _ = assert (plant_min_id <= plant_min_id_rhs)
let _ = assert (plant_crop_n = 24)
let _ = assert (worm_cells < worm_cells_rhs)
let _ = assert (ciona_cells < ciona_cells_rhs)
let _ = assert (platynereis_cells < platynereis_cells_rhs)
let _ = assert (male_neurons < male_neurons_rhs)
let _ = assert (banc_neurons < banc_neurons_rhs)
let _ = assert (ciona_gaba_flag = 0)
let _ = assert (platynereis_gaba_flag = 0)
