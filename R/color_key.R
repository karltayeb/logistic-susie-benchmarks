# colors for each method, keep consistent across
plot_colors <- function(){
  c(
    # SERs 
    'glm_ser' = 'aquamarine1',
    'glm_ser2' = 'aquamarine3',
    'uvb_ser' = 'brown1',
    'vb_ser' = 'darkgoldenrod1',
    'linear_ser' = 'royalblue4',
    'poly_ser_M2' = '#088293',
    #'uvb_ser2' = 'brown3', 
    #'glm_ser_aug' = '#088465', 
    #'vb_ser2' = 'darkgoldenrod3',
    #'vb_ser_corrected' = 'darkorange1',
    #'veb_ser' = 'darkorchid1',
    #'quad_ser' = 'deeppink1',
    #'uvb_ser_re' = 'palevioletred3',
    #'bayes_ser' = '#5a5252',

    # SuSiE
    'linear_susie_L5' =  'royalblue4',
    'binsusie_L5' = 'olivedrab3',
    'ibss_uvb_L5' = 'brown1', 
    'poly_susie_M2_L5'= '#ea97f0',
    'poly_susie_M6_L5'= '#ea5ff4',
    'poly_susie_M10_L5'= '#e415f3',

    # 'ibss_glm_L5' = 'aquamarine1', 
    # 'ibss_glm2_L5' = 'aquamarine3', 
    # 'ibss_glm_aug_L5' = '#088465',
    # 'ibss_uvb2_L5' = 'brown3', 
    # 'ibss_vb_L5' = 'darkgoldenrod1', 
    # 'ibss_vb2_L5' = 'darkgoldenrod3',
    # 'ibss_vbc_L5' = 'darkorange1',
    # 'binsusie2_L5' = 'olivedrab4',
    # 'ibss2m_L5' = 'brown3',
    # 'ibss2m_L5_jax' = 'brown4',
    'none' = 'black'
  )
}
