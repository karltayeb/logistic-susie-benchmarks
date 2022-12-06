# Website ------
# define targets for building the website
# the notebooks for building the site are in notebooks directory
website <- list(
  tar_render(web_home, 'notebooks/index.Rmd', output_dir='docs'),
  tar_render(web_calibration, 'notebooks/calibration.Rmd', output_dir='docs')
)
