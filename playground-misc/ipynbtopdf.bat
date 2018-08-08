@echo off

for %%f in (*.ipynb) do jupyter nbconvert "%%f" --to=pdf --TemplateExporter.exclude_output=True
