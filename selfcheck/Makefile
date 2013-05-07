tests = test1 test_diaglist

.PHONY += selfcheck $(tests)

selfcheck : $(tests)
	@echo Done testing!

test1 : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test1/
test_diaglist : include_selfcheck:=$(FC_MODINC)$(selfcheck_dir)test_diaglist/
test% : clean model
	$(SHELL) $(selfcheck_dir)runselfcheck.sh $@
