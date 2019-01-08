from ij.gui import GenericDialog, DialogListener
from java.awt import CheckboxGroup, Checkbox, Panel

def set_numeric(val):
	print("Numeric value is " + str(val));

def set_string(val):
	print("String value is " + str(val));

def set_choice(val):
	print("Choice is: "  + str(val));

def set_radio_button(val):
	print("Radio button is set to " + str(val));

def set_checkbox(val):
	print("Checkbox value is set to " + str(val));

def set_radiobutton(val):
	print("Radiobuttongroup is set to "  + str(val));
	
class MyControlDefinition:
	Numeric, String, Checkbox, Choice, RadioButtonGroup = range(5);
	def __init__(self, label, control_type, default_val, underlying_data_setter, choices=None, enabled=True):
		self.label = str(label);
		self.control_type = int(control_type);
		if self.control_type < 0 or self.control_type > MyControlDefinition.RadioButtonGroup:
			raise NotImplementedError();
		self.default_val = default_val;
		self.setter = underlying_data_setter;
		self.enabled = enabled;
		self.choices = choices;
		self.checkboxes = [];

	def addControl(self, dialog):
		if self.control_type==MyControlDefinition.Numeric:
			dialog.addNumericField(self.label, self.default_val, 2);
			dialog.getNumericFields()[-1].setEnabled(self.enabled);
		elif self.control_type==MyControlDefinition.String:
			dialog.addStringField(self.label, self.default_val);
			dialog.getStringFields()[-1].setEnabled(self.enabled);
		elif self.control_type==MyControlDefinition.Choice:
			if not self.default_val in self.choices:
				raise IndexError("default value isn''t in list of choices");
			dialog.addChoice(self.label, self.choices, self.default_val);
			dialog.getChoices()[-1].setEnabled(self.enabled);
		elif self.control_type==MyControlDefinition.Checkbox:
			if not isinstance(self.default_val, bool):
				raise TypeError();
			dialog.addCheckbox(self.label, self.default_val);
			dialog.getCheckboxes()[-1].setEnabled(self.enabled);
		elif self.control_type==MyControlDefinition.RadioButtonGroup:
			if not self.default_val in self.choices:
				raise IndexError("default value isn''t in list of choices");
			panel = Panel();
			cbg = CheckboxGroup();
			for ch in self.choices:
				cb = Checkbox(ch, cbg, ch==self.default_val)
				cb.setEnabled(self.enabled);
				self.checkboxes.append(cb);
				panel.add(cb);
			dialog.addPanel(panel);


	def setValue(self, value):
		self.setter(value);

controls = []

controls.append(MyControlDefinition("numeric", MyControlDefinition.Numeric, 3.3, set_numeric));
controls.append(MyControlDefinition("string", MyControlDefinition.String, "hello", set_string, enabled=False));
controls.append(MyControlDefinition("another numeric", MyControlDefinition.Numeric, 1, set_numeric));
controls.append(MyControlDefinition("a choice: ", MyControlDefinition.Choice, "choice 1", set_choice, choices=["choice 1", "choice 2", "choice 3"]));
controls.append(MyControlDefinition("a checkbox", MyControlDefinition.Checkbox, True, set_checkbox));
controls.append(MyControlDefinition("a radiobuttongroup", MyControlDefinition.RadioButtonGroup, "three", set_radiobutton, choices=["one", "two", "three"], enabled=False));

dialog = GenericDialog("Test");
for control in controls:
	control.addControl(dialog);
dialog.showDialog();

print([c.label for c in controls]);

numeric_fields = dialog.getNumericFields()
print(numeric_fields);
numeric_controls = [c for c in controls if c.control_type==MyControlDefinition.Numeric];
for nc, nf in zip(numeric_controls, dialog.getNumericFields()):
	nc.setter(float(nf.getText()));
string_controls = [c for c in controls if c.control_type==MyControlDefinition.String];
for sc, sf in zip(string_controls, dialog.getStringFields()):
	sc.setter(sf.getText());
choice_controls = [c for c in controls if c.control_type==MyControlDefinition.Choice];
for cc, cf in zip(choice_controls, dialog.getChoices()):
	cc.setter(cf.getSelectedItem());
checkbox_controls = [c for c in controls if c.control_type==MyControlDefinition.Checkbox];
for cbc, cbf in zip(checkbox_controls, dialog.getCheckboxes()):
	cbc.setter(cbf.getState());
radiobuttongroup_controls = [c for c in controls if c.control_type==MyControlDefinition.RadioButtonGroup];
for rbc in radiobuttongroup_controls:
	rbc.setter(rbc.checkboxes[[cb.getState() for cb in rbc.checkboxes].index(True)].getLabel());
