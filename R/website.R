# Website ------
# define targets for building the website
# the notebooks for building the site are in notebooks directory
website <- list(
  tar_render(web_linear_fail, 'notebooks/linear_fail.Rmd', output_dir='docs'),
  tar_render(web_home, 'notebooks/index.Rmd', output_dir='docs'),
  tar_render(web_constant_sim, 'notebooks/constant_sim_results.Rmd', output_dir='docs'),
  tar_render(web_half_normal_sim, 'notebooks/half_normal_sim_results.Rmd', output_dir='docs'),
  tar_render(web_tilted_vs_jj, 'notebooks/tilted_vs_jj.Rmd', output_dir='docs'),
  tar_render(web_stochastic_susie, 'notebooks/stochastic_susie.Rmd', output_dir='docs'),
  tar_render(web_glm_abf_exploration, 'notebooks/glm_abf_exploration.Rmd', output_dir='docs'),
  tar_render(web_univariate_bayes_approximations, 'notebooks/univariate_bayes_lr_comparison.Rmd', output_dir='docs'),
  tar_render(web_yusha_ibss_pdac_factor3, 'notebooks/yusha_ibss_pdac_factor3.Rmd', output_dir='docs')
)
