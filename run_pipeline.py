"""run_pipeline.py - executa as etapas disponiveis em ordem."""
import importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
STEPS = ["01_download.py", "02_frequencies.py", "05_brazil.py",
         "04_report.py", "viz_all.py", "03_figures.py"]


def run_step(fn):
    spec = importlib.util.spec_from_file_location(fn[:-3], os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main()


def main():
    for fn in STEPS:
        print(f"\n===== {fn} =====")
        run_step(fn)
    print("\n>> Pipeline (etapas atuais) concluido.")


if __name__ == "__main__":
    main()
