from java.awt import CheckboxGroup, Checkbox, Panel

class MyControlDefinition:
	Numeric, String, Checkbox, Choice, RadioButtonGroup = range(5);
	def __init__(self, label, control_type, default_val, underlying_data_setter, choices=None, enabled=True):
		self.label = label;
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