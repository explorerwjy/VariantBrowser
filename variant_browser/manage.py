from flask.ext.script import Manager
from variant_browser import app
import variant_browser

manager = Manager(app)

@manager.command
def hello():
    print "hello"

@manager.command
def load_db():
    variant_browser.load_db()

@manager.command
def load_variants_file():
    variant_browser.load_variants_file()


@manager.command
def reload_variants():
    variant_browser.load_variants_file()

@manager.command
def load_gene_models():
    variant_browser.load_gene_models()


if __name__ == "__main__":
    manager.run()
