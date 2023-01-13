# Website ------
# define targets for building the website
# the notebooks for building the site are in notebooks directory
website <- list(
  tar_render(web_linear_fail, 'notebooks/linear_fail.Rmd', output_dir='docs'),
  tar_render(web_home, 'notebooks/index.Rmd', output_dir='docs'),
  tar_render(web_constant_sim, 'notebooks/constant_sim_results.Rmd', output_dir='docs'),
  tar_render(web_half_normal_sim, 'notebooks/half_normal_sim_results.Rmd', output_dir='docs'),
  tar_render(web_tilted_vs_jj, 'notebooks/tilted_vs_jj.Rmd', output_dir='docs')
)


library(dplyr)
tar_meta() %>%
  filter(stringr::str_detect(name, 'constant_fit')) %>%
  head()
