# v0.30.25
+ bump to get new release for CI/wheel fixes. thanks to @graphenn (#286) and @oyvinev (#297)

# v0.30.24
+ fix memory leak in raw_header (#281)

# v0.30.23
+ set samples to empty when given empty list.
+ Prevent memory leak if user calls VCF.close() (#280)

# v0.30.19
+ raise exception (instead of segfault) when calling num_called, etc. on variant without genotypes (#262 from @davmlaw)

# v0.30.15
+ fix num_unknown when gts012=True (#236)

# v0.30.14
+ use warnings instead of sys.stderr (#229 from @grahamgower)
+ use libdeflate in wheel build (#231 from @grahamgower)
+ use pytest instead of nose and update numpy stuff (#232 from @grahamgower)

# v0.30.13
+ fixes for mixed ploidy samples affecting `variant.num_het`,
  `variant.num_hom_ref`, etc. Also affects `variant.genotypes` and
  `variant.genotype` (see #227. Thanks @davmlaw and  @grahamgower for
  test-cases and insight.)

# v0.30.12
+ add variant.FILTERS (by @tomwhite, #149)

# v0.30.11
+ bump for CI

# v0.30.10
+ fix is_indel for <NON_REF> (see #217)

# 0.3.8
+ just bumping for CI+wheels

# 0.30.7
+ fix gt_types for phased variants (see: #198)
