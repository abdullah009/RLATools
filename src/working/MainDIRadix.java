package working;

public class MainDIRadix
{
	
	
	public static void main(String[] args)
	{
		System.out.println(args.length);
		if (args.length < 1) {
			System.out.println("Command line format: java -jar path/RLT_S.jar path/config_file_name.xml\n");
			System.exit(1);
		}
		//System.out.println("Record Linkage using single linkage clustering  " + args[0]);
		
		//System.out.println("Path: " + (new java.io.File("").getAbsolutePath()) + " Path: " + System.getProperty("user.dir"));
		
		initApp(args);
	}
	
	public static void initApp(String[] args)
	{
		DIRadix diRadix	= new DIRadix(args[0]); // args[0] => directory containing xml
		diRadix.clusterData();
	}
	
}
