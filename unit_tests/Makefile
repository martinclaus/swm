
unit_tests: testdriver

.PHONY += unit_tests

COMPILE_TEST_EXE = $(FC) $(FFLAGS) $(CPPFLAGS) $(FC_MODINC)$(builddir) $(FC_MODINC)$(includedir) $(FC_MODOUT)$(builddir) -o$(bindir)$@ $< $(patsubst %,$(builddir)%,$(filter %.o,$^))

testdriver$(EXEEXT) : $(utest_dir)test_driver.f90 utest_list.$(OBJEXT) utest_app.$(OBJEXT) utest_io.$(OBJEXT) list_mod.$(OBJEXT) app.$(OBJEXT) io_module.$(OBJEXT) testing.$(OBJEXT) calendar_module.$(OBJEXT) f_udunits_2.$(OBJEXT) logging.$(OBJEXT) component_module.$(OBJEXT)
	$(COMPILE_TEST_EXE) $(LDFLAGS)
	$(bindir)$@

testing.o : $(utest_dir)testing.f90
	@$(COMPILE)

utest_list.$(OBJEXT) : $(utest_dir)utest_list.f90 testing.$(OBJEXT) list_mod.$(OBJEXT)
	@$(COMPILE)

utest_app.$(OBJEXT) : $(utest_dir)utest_app.f90 testing.$(OBJEXT) app.$(OBJEXT) component_module.$(OBJEXT)
	@$(COMPILE)

utest_io.$(OBJEXT) : $(utest_dir)utest_io.f90 testing.$(OBJEXT) io_module.$(OBJEXT) list_mod.$(OBJEXT)
	@$(COMPILE)
