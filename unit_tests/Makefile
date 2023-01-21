
unit_tests: testdriver

.PHONY += unit_tests

COMPILE_TEST_EXE = $(FC) $(FFLAGS) $(CPPFLAGS) $(FC_MODINC)$(builddir) $(FC_MODINC)$(includedir) $(FC_MODOUT)$(builddir) -o$(bindir)$@ $< $(patsubst %,$(builddir)%,$(filter %.o,$^))

testdriver$(EXEEXT) : $(utest_dir)test_driver.f90 utest_list.$(OBJEXT) utest_app.$(OBJEXT) list_mod.$(OBJEXT) app.$(OBJEXT) testing.$(OBJEXT)
	@$(COMPILE_TEST_EXE)
	@$(bindir)$@

testing.o : $(utest_dir)testing.f90
	@$(COMPILE)

utest_list.$(OBJEXT) : $(utest_dir)utest_list.f90 testing.$(OBJEXT) list_mod.$(OBJEXT)
	@$(COMPILE)

utest_app.$(OBJEXT) : $(utest_dir)utest_app.f90 testing.$(OBJEXT) app.$(OBJEXT)
	$(COMPILE)
