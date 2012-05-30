rm outputs/*zjetsall.root;
hadd -f outputs/PUsysUp_0_zjetsall.root outputs/PUsysUp_0_zjets*;
hadd -f outputs/PUsysUp_1_zjetsall.root outputs/PUsysUp_1_zjets*;
hadd -f outputs/PUsysUp_2_zjetsall.root outputs/PUsysUp_2_zjets*;

hadd -f outputs/PUsysUp_0_twdr.root outputs/PUsysUp_0_tw_dr.root outputs/PUsysUp_0_atw_dr.root;
hadd -f outputs/PUsysUp_1_twdr.root outputs/PUsysUp_1_tw_dr.root outputs/PUsysUp_1_atw_dr.root;
hadd -f outputs/PUsysUp_2_twdr.root outputs/PUsysUp_2_tw_dr.root outputs/PUsysUp_2_atw_dr.root;

hadd -f outputs/PUsysUp_0_others.root outputs/PUsysUp_0_ww.root outputs/PUsysUp_0_wz.root outputs/PUsysUp_0_zz.root outputs/PUsysUp_0_t.root outputs/PUsysUp_0_at.root outputs/PUsysUp_0_ts.root outputs/PUsysUp_0_ats.root outputs/PUsysUp_0_wjets.root outputs/PUsysUp_0_qcd_mu.root ;
hadd -f outputs/PUsysUp_1_others.root outputs/PUsysUp_1_ww.root outputs/PUsysUp_1_wz.root outputs/PUsysUp_1_zz.root outputs/PUsysUp_1_t.root outputs/PUsysUp_1_at.root outputs/PUsysUp_1_ts.root outputs/PUsysUp_1_ats.root outputs/PUsysUp_1_wjets.root outputs/PUsysUp_1_qcd_mu.root ;
hadd -f outputs/PUsysUp_2_others.root outputs/PUsysUp_2_ww.root outputs/PUsysUp_2_wz.root outputs/PUsysUp_2_zz.root outputs/PUsysUp_2_t.root outputs/PUsysUp_2_at.root outputs/PUsysUp_2_ts.root outputs/PUsysUp_2_ats.root outputs/PUsysUp_2_wjets.root outputs/PUsysUp_2_qcd_mu.root ;

hadd -f outputs/PUsysUp_0_others_2.root outputs/PUsysUp_0_others.root outputs/PUsysUp_0_zjetsall.root
hadd -f outputs/PUsysUp_1_others_2.root outputs/PUsysUp_1_others.root outputs/PUsysUp_1_zjetsall.root
hadd -f outputs/PUsysUp_2_others_2.root outputs/PUsysUp_2_others.root outputs/PUsysUp_2_zjetsall.root
