// 
// This defines the abstract class Process from which every user-defined process
// will inherit.
public abstract class Process {
	// Constructor
	protected Process(movieData owner, String name, String dateTime)
	{
		owner_ = owner;
	}
	
	// Get the name of the process
	public String getName() {
		return name_;
	}
	
	// Get the date when the process has been finished
	public String getDateTime() {
		return dateTime_;
	}		
	
	// Get the comment the user could have wrote on the process GUI.
	public String getComment() {
		return comment_;
	}
	
	// abstract methods section
	
	// Get the help text
	public abstract String getHelp();
	
	// make a sanity check of the process
	public abstract void sanityCheck(boolean full);
	
	// more abstract methods...
	
	// protected field section
	protected const movieData owner_;
	protected String dateTime_;	
	protected String comment_;
	
	// private fields
	private String name_;
}
