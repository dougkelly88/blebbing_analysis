from ij.gui import GenericDialog, DialogListener

def set_hello(val):
	print("Hello entry is " + str(val));

def set_another(val):
	print("Another numeric is " + str(val));

class MyControlDefinition:
	Numeric, String, Checkbox, Choice, RadioButtonGroup = range(5);
	def __init__(self, label, control_type, default_val, setter, choices=None, enabled=True):
		self.label = str(label);
		self.control_type = int(control_type);
		if self.control_type < 0 or self.control_type > MyControlDefinition.RadioButtonGroup:
			raise NotImplementedError();
		self.default_val = default_val;
		self.setter = setter;
		self.enabled = enabled;

	def addControl(self, dialog):
		if self.control_type==MyControlDefinition.Numeric:
			dialog.addNumericField(self.label, self.default_val, 2);
			dialog.getNumericFields()[-1].setEnabled(self.enabled);
		elif self.control_type==MyControlDefinition.String:
			dialog.addStringField(self.label, self.default_val);
			dialog.getStringFields()[-1].setEnabled(self.enabled);

	def setValue(self, value):
		self.setter(value);

controls = []

controls.append(MyControlDefinition("hello", MyControlDefinition.Numeric, 3.3, set_hello));
controls.append(MyControlDefinition("goodbye", MyControlDefinition.String, "forever?", set_another, enabled=False));
controls.append(MyControlDefinition("another numeric", MyControlDefinition.Numeric, 1, set_another));

dialog = GenericDialog("Test");
for control in controls:
	control.addControl(dialog);
dialog.showDialog();

print([c.label for c in controls]);

numeric_fields = dialog.getNumericFields()
print(numeric_fields);
numeric_controls = [c for c in controls if c.control_type==MyControlDefinition.Numeric];
for nc, nf in zip(numeric_controls, numeric_fields):
	nc.setter(float(nf.getText()));
	