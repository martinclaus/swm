tests = utest_str test1 test_diaglist

.PHONY += selfcheck $(tests)

selfcheck : $(tests)
	@echo Done testing!

test1 : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test1/
test_diaglist : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_diaglist/
test_nonlinear : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_nonlinear/
test_nl_tracer : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_nl_tracer/
test_oscfrc : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_oscfrc/
test% : clean model
	@echo $@
	$(SHELL) $(selfcheck_dir)runselfcheck.sh $@

utest_str : $(selfcheck_dir)utest_str/utest_str.f90 str.o
	$(FC) $(FFLAGS) $(CPPFLAGS) $(FC_MODINC)$(builddir) $(builddir)str.o $< -o $(selfcheck_dir)$@/$@
	$(selfcheck_dir)$@/$@
