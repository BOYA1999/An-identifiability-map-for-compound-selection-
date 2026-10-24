PY ?= python
PKG := package/ExpansionRx_reproducibility

.PHONY: verify replay derived checksums clean

verify:
	cd $(PKG) && $(PY) -B -X utf8 verify_package.py

replay:
	cd $(PKG) && $(PY) -B -X utf8 replay_saved.py

derived:
	cd $(PKG) && $(PY) -B -X utf8 reproduce_derived.py

checksums:
	cd $(PKG) && $(PY) -B -X utf8 -c "import json,hashlib,pathlib; m=json.load(open('MATERIALS_MANIFEST.json')); print('\n'.join(hashlib.sha256((pathlib.Path('.')/e['path']).read_bytes()).hexdigest()+'  '+e['path'] for e in m['files']))" > ../../checksums/SHA256SUMS

clean:
	rm -rf $(PKG)/replay_output $(PKG)/regeneration_output
