spec_dec_dict = {
"01_Z_m_p"   :  "B0 -> { D- -> K+ pi- pi- } { D+ -> K- pi+ pi+ } {K*0 -> K+ pi-}",
"02_Z_m_p"   :  "B0 -> { D*- -> {D- -> K+ pi- pi-} pi0 } { D+ -> K- pi+ pi+ } { K*0 -> K+ pi-}",
"02_P_z_p"   :  "B0 -> { D*- -> {D0b -> K+ pi- } pi- } { D+ -> K- pi+ pi+ } { K*0 -> K+ pi-}",
"04_Z_m_p"   :  "B0 -> { D*- -> {D- -> K+ pi- pi-} pi0 } { D*+ -> {D+ -> K- pi+ pi+} pi0 } { K*0 -> K+ pi-}",
"04_P_z_p"   :  "B0 -> { D*- -> {D0b -> K+ pi-} pi- } { D*+ -> {D+ -> K- pi+ pi+} pi0 } { K*0 -> K+ pi-}",
"04_Z_z_z"   :  "B0 -> { D*- -> {D0b -> K+ pi-} pi- } { D*+ -> {D0 -> K- pi+} pi+ } { K*0 -> K+ pi-}",
"04_P_z_pst" :   "B0 -> { D*- -> {D0b -> K+ pi-} pi- } { D*+ -> {D0 -> K- pi+} pi+ } { K*0 -> K+ pi-}",
"05_P_z_p"   :   "B+ -> { D0b -> K+ pi- } { D+ -> K- pi+ pi+ } {K*0 -> K+ pi-}",
"06_P_z_p"   :   "B+ -> { D*0b -> {D0b -> K+ pi-} pi0 } { D+ -> K- pi+ pi+ } {K*0 -> K+ pi-}",
"07_P_z_p"    :  "B+ -> { D0b -> K+ pi- } { D*+ -> {D+ -> K- pi+ pi+} pi0 } {K*0 -> K+ pi-}",
"07_Z_z_z"    :  "B+ -> { D0b -> K+ pi- } { D*+ -> {D0 -> K- pi+} pi+ } {K*0 -> K+ pi-}",
"07_P_z_pst"   : "B+ -> { D0b -> K+ pi- } { D*+ -> {D0 -> K- pi+} pi+ } {K*0 -> K+ pi-}",
"08_P_z_p"     : "B+ -> { D*0b -> {D0b -> K+ pi-} pi0 } { D*+ -> {D+ -> K- pi+ pi+} pi0 } {K*0 -> K+ pi-}",
"08_Z_z_z"     : "B+ -> { D*0b -> {D0b -> K+ pi-} pi0 } { D*+ -> {D0 -> K- pi+} pi+ } {K*0 -> K+ pi-}",
"08_P_z_pst"   : "B+ -> { D*0b -> {D0b -> K+ pi-} pi0 } { D*+ -> {D0 -> K- pi+} pi+ } {K*0 -> K+ pi-}",
"09_Z_z_z"     : "B0 -> { D0b -> K+ pi- } { D0 -> K- pi+ } {K*0 -> K+ pi-}",
"10_Z_z_z"    :  "B0 -> { D*0b -> {D0b -> K+ pi-} pi0 } { D0 -> K- pi+ } {K*0 -> K+ pi-}",
"12_Z_z_z"    :  "B0 -> { D*0b -> {D0b -> K+ pi-} pi0 } { D*0 -> {D0 -> K- pi+} pi0 } {K*0 -> K+ pi-}",
"norm7"       :  "B0 -> { D0b -> K+ pi- } { D0 -> K- pi+ pi+ pi- } K+",
"norm8"       :  "B0 -> { D- -> K+ pi- pi- } { D0 -> K- pi+ pi+ pi- } K+",}

for key, value in spec_dec_dict.items():
    file = open(f"decay_files/{key}.decay", "w")
    file.write(value)
    file.close()
