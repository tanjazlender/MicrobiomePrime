import configparser

class Config:
    def __init__(self, config_file='scripts/variables.ini'):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)
        
        if 'settings' not in self.config:
            raise KeyError("The 'settings' section is missing in the configuration file.")
        
        settings = self.config['settings']
        
        def safe_get(key):
            value = settings.get(key, None)
            if value is None:
                raise KeyError(f"Missing required configuration key: '{key}'")
            return value

        def convert(value):
            value = value.strip()
            if value == "":
                return None  # Handle empty values if needed
            try:
                # Convert to int if possible
                int_value = int(value)
                return int_value
            except ValueError:
                # Convert to float if not int
                try:
                    float_value = float(value)
                    return int(float_value) if float_value.is_integer() else float_value
                except ValueError:
                    # Return as a string if conversion fails
                    return value

        # Extract configuration settings
        self.kmer_sensitivity_cutoff = convert(safe_get('kmer_sensitivity_cutoff'))
        self.kmer_specificity_cutoff = convert(safe_get('kmer_specificity_cutoff'))
        self.marker_sensitivity_cutoff = convert(safe_get('marker_sensitivity_cutoff'))
        self.marker_specificity_cutoff = convert(safe_get('marker_specificity_cutoff'))
        self.min_amplicon_length = convert(safe_get('min_amplicon_length'))
        
        # Ensure `raw_target` is a string
        self.raw_target = settings.get('target', '').strip()
        self.target_group_name = settings.get('target_group_name', '')
        self.specificity_exception = [e.strip() for e in settings.get('specificity_exception', '').split(',') if e.strip()]

        # Process target information
        self.target = [x.strip() for x in self.raw_target.split(',') if x.strip()]
        self.target_combined = " ".join(self.target)
        self.target_group_ID = self.target_group_name if self.target_group_name else self.target_combined.replace(" ", "-")

    def __repr__(self):
        return (f"Config(kmer_sensitivity_cutoff={self.kmer_sensitivity_cutoff}, "
                f"kmer_specificity_cutoff={self.kmer_specificity_cutoff}, "
                f"marker_sensitivity_cutoff={self.marker_sensitivity_cutoff}, "
                f"marker_specificity_cutoff={self.marker_specificity_cutoff}, "
                f"min_amplicon_length={self.min_amplicon_length}, "
                f"raw_target={self.raw_target}, "
                f"target={self.target}, "
                f"target_combined={self.target_combined}, "
                f"target_group_name={self.target_group_name}, "
                f"target_group_ID={self.target_group_ID}, "
                f"specificity_exception={self.specificity_exception})")