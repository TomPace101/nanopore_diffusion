
import simproc.requesthandler.debug as debug
import simproc.requesthandler.yaml_manager as yaml_manager

class AnotherDummyShellRequest(debug.DummyShellRequest):
  pass

yaml_manager.register_classes([AnotherDummyShellRequest])