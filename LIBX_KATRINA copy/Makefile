all:
	./1_create_readme_proto.py
	./2_get_unique_decorations.py
	./3_get_oxidations.py
	run_multish --np=6 --file=commands_decorations.sh
	bash commands_decorations.sh
	-run_multish --np=12 --file=commands_oxidations.sh
	-bash commands_oxidations.sh
	./4_grab_protos.py
	mv commands_aflowin.sh commands_aflowin_KEEP.sh
	./4_grab_protos.py --debug 2>&1 | tee out_4_grab_protos_`date +%Y%m%d`.txt
	xz out_4_grab_protos_`date +%Y%m%d`.txt
	cmp --silent commands_aflowin_KEEP.sh commands_aflowin.sh

