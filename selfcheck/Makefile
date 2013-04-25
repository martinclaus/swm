tests = test1 test_diaglist

.PHONY += selfcheck $(tests)

selfcheck : $(tests)
	@echo Done testing!

test1 : include_modelHeader:=$(selfcheck_dir)test1/
test_diaglist : include_modelHeader:=$(selfcheck_dir)test_diaglist/
test% : model
	$(SHELL) $(selfcheck_dir)runselfcheck.sh $@
