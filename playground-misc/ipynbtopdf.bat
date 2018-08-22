@echo off

for %%f in (*.ipynb) do (
	echo ========BEGIN: Convert %%f
	jupyter nbconvert "%%f" --to=pdf --TemplateExporter.exclude_output=True
	echo ========FINISH: Convert %%f
	echo 
)