tests = utest_str test1 test_diaglist

.PHONY += selfcheck $(tests)

selfcheck : $(tests)
	@echo Done testing!

selftest_1 : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_1/
selftest_diaglist : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_diaglist/
selftest_nonlinear : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_nonlinear/
selftest_nl_tracer : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_nl_tracer/
selftest_oscfrc : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_oscfrc/
selftest_performance : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_performance/
selftest_logging : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_logging/
selftest_fixed_length_output : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_fixed_length_output/

selftest_% : clean model
	@echo run $@ with $(bindir)model
	$(SHELL) $(selfcheck_dir)runselfcheck.sh $@ $(bindir)model

utest_str : $(selfcheck_dir)utest_str/utest_str.f90 str.o
	$(FC) $(FFLAGS) $(CPPFLAGS) $(FC_MODINC)$(builddir) $(builddir)str.o $< -o $(selfcheck_dir)$@/$@
	$(selfcheck_dir)$@/$@
