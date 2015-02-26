/*
 * Abdullah-Al Mamun
 */

package working;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;


public class DIRadix
{
	int lmerTotal		= 4;
	int lmerSize		= 3; // according to config
	int alphabetType	= 0; // according to config
	final int ALPHABET_SIZE		= 256;
	final int[] ALPHABET_SIZE_LIST	= {26, 10, 36, 256};
	final int RECORD_TOTAL		= 500000;
	
	SimpleDateFormat dateFormat;
	long readT, writeT, sortT, clusterExactT, blockT, edgeListT, clusterApproxT, sccT, clusterFinalT, totalT, startT;
	int threshold, pointTotal, lenMax, clusterTotal;
	String configFileStr, outFileAddrStr, inFileAddrStr, strSample, fieldStr;
		
	private ArrayList<ArrayList<String>> pointArr;
	private ArrayList<ArrayList<Integer>> clusterExactIndArr, indexDatasetArr;
	private ArrayList<ArrayList<ArrayList<String>>> clusterArr;
	private ArrayList<int[]> edgeArr, editAtt, truncAtt;
	private ArrayList<ArrayList<int[]>> reversalAtt;
	private ArrayList<String> fileNameArr;
	private ArrayList<Integer> truncRate, compArr, fieldTotalArr, attrCommonArr;
	private StrForSort[] strDataArr;
	
	
	public DIRadix(String configFileStr)
    {
		this.configFileStr	= configFileStr;
		initThis();
    }
	
	private void initThis()
	{
		indexDatasetArr	= new ArrayList<ArrayList<Integer>>();
		fieldTotalArr	= new ArrayList<Integer>();
		compArr 		= new ArrayList<Integer>();
		editAtt 		= new ArrayList<int[]>();
		truncAtt 		= new ArrayList<int[]>();
		truncRate		= new ArrayList<Integer>();
		reversalAtt 	= new ArrayList<ArrayList<int[]>>();
		
		fieldStr		= getTitle();
		fileNameArr		= getInputFileNameList();
		attrCommonArr	= getCommonAttrList();
		threshold		= getInputThreshold();
		getInputComparisonPara();

		pointArr		= new ArrayList<ArrayList<String>>();
						
		//System.out.println(threshold + " : " + fileNameArr.get(0) + " " + configFileStr);
		dateFormat		= new SimpleDateFormat("hh:mm:ss");
		long currTS1	= System.currentTimeMillis();
		readDataFromFile();
		long currTS2	= System.currentTimeMillis();
		readT			= currTS2 - currTS1;
		//System.out.println("readFromFile: " + readT + " ms");
		
		/*
		for(int i = 0; i < editAtt.size(); ++i)
		{
			for(int j = 0; j < editAtt.get(i).length; ++j)
				System.out.print(editAtt.get(i)[j] + "\t");
			System.out.println(":" + i);
		}
		*/
	}
	
	private void readDataFromFile()
	{
		int i, j;//, pointPerFile, count;
		boolean doneF;
				
		for(int l = 0; l < fileNameArr.size(); ++l)
		{
	        try
	        {
	        	BufferedReader reader;
	        	ArrayList<String> rowStrArr;
	        	String lineStr, tempStr;
	        	String[] attrArr;
	        		        
	        	lineStr		= null;
	        	reader 		= new BufferedReader(new FileReader(fileNameArr.get(l)));
	        	doneF		= false;
	        	
	        	//lineStr = reader.readLine(); // for field row
		        while((lineStr = reader.readLine()) != null)
		        {
		        	if(lineStr.length() > 10)
		        	{	
			        	rowStrArr	= new ArrayList<String>();
			        	tempStr		= lineStr;
			        	attrArr		= tempStr.split(",");
			        	for(i = 0; i < attrArr.length; ++i)
			        		rowStrArr.add(attrArr[i].toLowerCase().trim());
			        	rowStrArr.add(fileNameArr.get(l));
			        	rowStrArr.add(Integer.toString(l));
			        	pointArr.add(rowStrArr);  
			        	if(!doneF)
			        	{
			        		fieldTotalArr.add(attrArr.length);
			        		doneF	= true;
			        	}
			        	//++count;
			        	//if(count >= pointPerFile) break;
		        	}
		        }
		        
		        reader.close();
		       
	        }
	        catch(FileNotFoundException e)
	        {
		        // TODO Auto-generated catch block
		        e.printStackTrace();
	        }
	        catch(IOException e)
	        {
	        	e.printStackTrace();
	        }
		}
        
		pointTotal	= pointArr.size();
        System.out.println("total record: " + pointTotal);
        
        String strCC;
        int lenDiff, attrCommonTotal;
        
        strSample	= "00000000000000000000000000000000000000000000";
        lenMax		= 0;
        strDataArr	= new StrForSort[pointTotal];
        
        attrCommonTotal	= attrCommonArr.size();
        
        for(i = 0; i < pointTotal; ++i)
        {
        	strCC	= "";
        	for(j = 0; j < attrCommonTotal; ++j)
        	{
        		strCC		+= pointArr.get(i).get(indexDatasetArr.get(Integer.parseInt(pointArr.get(i).get(pointArr.get(i).size() - 1))).get(attrCommonArr.get(j)));
        		if(j < (attrCommonTotal - 1))
        			strCC	+= "#";
        	}
        	strDataArr[i]	= new StrForSort(i, strCC);
        	if(strCC.length() > lenMax)
        		lenMax		= strCC.length();
        	//System.out.println(i + ":" + strCC + " len: " + strCC.length());
        }
        for(i = 0; i < pointTotal; ++i)
        {
        	lenDiff			= lenMax - strDataArr[i].str.length();
        	if(lenDiff > 0)
        		strDataArr[i].str	= strDataArr[i].str.concat(strSample.substring(0, lenDiff));
        }
		
	}


	/*
	 * cluster records finding exact match, then approximate match and finally combining them  
	 */
	public void clusterData()
	{
		System.out.println("Clustering data");
		
		long currTS1, currTS2, currTS3;
		
		currTS1			= System.currentTimeMillis();
		startT			= currTS1;
		findExactCluster(); // find clusters using exact match
		clusterExactT	= System.currentTimeMillis() - currTS1;
		
		currTS2			= System.currentTimeMillis();
		ArrayList<Integer> rootArr	= findApproxCluster(); // find clusters using approximate match
		clusterApproxT	= System.currentTimeMillis() - currTS2;
		
		currTS3			= System.currentTimeMillis();
		findFinalCluster(rootArr); // combine exact and approximate match
		clusterFinalT	= System.currentTimeMillis() - currTS3;
		
		output(); // write output to a file
		
		totalT			= System.currentTimeMillis() - currTS1;
		
		System.out.println("Total Cluster " + clusterTotal);
		//System.out.println("Total Cluster " + clusterTotal + "\n" + "clusterExactT: " + clusterExactT + "\tclusterApproxT: " + clusterApproxT + "\tclusterFinalT: " + clusterFinalT + "\ttotalT: " + totalT);
	}
	
	
	/*
	 * find clusters with no errors by grouping records using exact match 
	 */
	private void findExactCluster()
	{
		sortData(); // sort records using some attributes as values (radix sort)
		
		int i;
		ArrayList<Integer> clusterRowArr;
				
		clusterExactIndArr	= new ArrayList<ArrayList<Integer>>();
		clusterRowArr		= new ArrayList<Integer>();
		clusterRowArr.add(strDataArr[0].ind);
		
		for(i = 1; i < pointTotal; ++i)
		{
			if(strDataArr[i].str.equalsIgnoreCase(strDataArr[i - 1].str))
				clusterRowArr.add(strDataArr[i].ind);
			else
			{
				clusterExactIndArr.add(clusterRowArr);
				clusterRowArr	= new ArrayList<Integer>();
				clusterRowArr.add(strDataArr[i].ind);
			}
		}
		clusterExactIndArr.add(clusterRowArr);
		
		//System.out.println("total exact clusters: " + clusterExactIndArr.size());
	}
	
	
	/*
	 * sort(radix sort) records using some attributes
	 */
	private void sortData()
	{
		int i, j, k;
		StrForSort[] tempArr;
		int[] countArr;
		
		tempArr	= new StrForSort[pointTotal];
		
		for(i = lenMax - 1; i >= 0; --i)
		{
			countArr	= new int[ALPHABET_SIZE];
			
			for(j = 0; j < pointTotal; ++j)
				countArr[strDataArr[j].str.charAt(i) + 1]++;
			
			for(k = 1; k < 256; ++k)
				countArr[k]	+= countArr[k - 1];
			
			for(j = 0; j < pointTotal; ++j)
				tempArr[countArr[strDataArr[j].str.charAt(i)]++]	= strDataArr[j];
			
			for(j = 0; j < pointTotal; ++j)
				strDataArr[j]	= tempArr[j];
		}
		
		tempArr		= null;
		countArr	= null;
	}
	

	/*
	 * find approximate cluster using approximate match by blocking, generating edgelist followed by findig connected components 
	 */
	private ArrayList<Integer> findApproxCluster()
	{
				
		createClusterEdgeList(createBlock()); // create blocks followed by generating edgelist within blocks 
		
		return findConnComp(); // find connected components on the graph generated by edgelist as edges and records as points
	}
	
	
	/*
	 * create blocks of records using LMER(here LMER = 3) characters of last name
	 */
	private ArrayList<ArrayList<Integer>> createBlock()
	{
		int i, j, k, blockTotal, strLen, blockID, indFieldDataset, indBlockField;
		ArrayList<ArrayList<Integer>> blockArr;
		ArrayList<Integer> tempArr;
		ArrayList<String> record;
		int[] coderecordArr;
		String blockFieldStr;
		
		indBlockField	= getBlockInd();
		String typeTest	= pointArr.get(0).get(indexDatasetArr.get(0).get(indBlockField));
		if (DIRadix.isNumeric(typeTest))
			alphabetType	= 1;
		
		if(alphabetType == 0)
			strSample	= "aaaaaaaaaaa";
		else
			strSample	= "0000000000";
		
		blockTotal	= (int) Math.pow(ALPHABET_SIZE_LIST[alphabetType], lmerSize);
		blockArr	= new ArrayList<ArrayList<Integer>>(blockTotal);
		//System.out.println(indBlockField + "\t" + alphabetType + "\t" + lmerSize + "\t" + threshold + "\t block total: " + blockTotal);
		for(i = 0; i < blockTotal; ++i)
		{
			tempArr	= new ArrayList<Integer>();
			blockArr.add(tempArr);
		}
		
		for(i = 0; i < clusterExactIndArr.size(); ++i)
		{
			record			= pointArr.get(clusterExactIndArr.get(i).get(0));
			indFieldDataset	= indexDatasetArr.get(Integer.parseInt(record.get(record.size() - 1))).get(indBlockField);
			if(indFieldDataset < 0)
				continue;
			blockFieldStr	= record.get(indFieldDataset);
			strLen			= blockFieldStr.length();
			if(strLen < lmerSize)
				blockFieldStr	= strSample.substring(0, lmerSize - strLen).concat(blockFieldStr);
			strLen			= blockFieldStr.length(); 
			coderecordArr	= new int[strLen];
			for(j = 0; j < blockFieldStr.length(); ++j)
			{
				if( ( (alphabetType == 0) && !((blockFieldStr.codePointAt(j) >= 97 && blockFieldStr.codePointAt(j) <= 122)) )
				|| ( (alphabetType == 1) && !((blockFieldStr.codePointAt(j) >= 48 && blockFieldStr.codePointAt(j) <= 57)) )
				|| ( (alphabetType == 2) && (!(((blockFieldStr.codePointAt(j) >= 97 && blockFieldStr.codePointAt(j) <= 122)) || ((blockFieldStr.codePointAt(j) >= 48 && blockFieldStr.codePointAt(j) <= 57)))) )
				)
				{
					blockFieldStr	= blockFieldStr.substring(0, j) + blockFieldStr.substring(j + 1);
					--j;
				}
				else
					coderecordArr[j]	= blockFieldStr.codePointAt(j);
			}
			
			//blockID=0;
			for(j = 0; j < blockFieldStr.length() - lmerSize + 1; ++j)
			{
				blockID	= 0;
				for(k = 0; k < lmerSize; ++k)
					if(alphabetType == 0)
							blockID	+= (coderecordArr[j + k] - 97) * (int) Math.pow(ALPHABET_SIZE_LIST[alphabetType], lmerSize - k - 1);
					else if(alphabetType == 1)
						blockID	+= (coderecordArr[j + k] - 48) * (int) Math.pow(ALPHABET_SIZE_LIST[alphabetType], lmerSize - k - 1);
					else if(alphabetType == 2)
					{
						if(coderecordArr[j + k] >= 97)
							blockID	+= (coderecordArr[j + k] - 97) * (int) Math.pow(ALPHABET_SIZE_LIST[alphabetType], lmerSize - k - 1);
						else
							blockID	+= (coderecordArr[j + k] - 22) * (int) Math.pow(ALPHABET_SIZE_LIST[alphabetType], lmerSize - k - 1); // 48 - 26
					}
				//System.out.println(i + ":" + blockFieldStr + " id:" + blockID + " len:" + blockFieldStr.length());
				if(!(blockID < 0 || blockID >= blockTotal))
					blockArr.get(blockID).add(i);
			}
		}
		
		return blockArr;
	
	}
	
	
	/*
	 * creates edge list of records within blocks
	 */
	private void createClusterEdgeList(ArrayList<ArrayList<Integer>> blockArr)
	{
		int i, blockTotal;
		
		edgeArr		= new ArrayList<int[]>();
		blockTotal	= blockArr.size();
		
		for(i = 0; i < blockTotal; ++i)
		{
			if(blockArr.get(i).size() > 0)
				generateEdgilist(blockArr.get(i));
		}
	}
	
	/*
	 * generate edge list within a block
	 */
	private void generateEdgilist(ArrayList<Integer> blockRowArr)
	{
		int i, j, n, temp, blockItemTotal;
		int[] vectorArr;
		ArrayList<Integer> tempArr;
		ArrayList<ArrayList<String>> dataArr;
		
		n				= 0;
		blockItemTotal	= blockRowArr.size();
		vectorArr		= new int[blockItemTotal];
		tempArr			= new ArrayList<Integer>();
		dataArr			= new ArrayList<ArrayList<String>>();
		
		for(i = 0; i < blockItemTotal; ++i)
			dataArr.add(pointArr.get(clusterExactIndArr.get(blockRowArr.get(i)).get(0)));
		
		for (i = 0; i < blockItemTotal; i++)
		{
			if (vectorArr[i] == 0)
			{
				++n;
				tempArr.add(i);
				vectorArr[i]	= n;
				while(tempArr.size()>0)
				{
					temp 		= tempArr.remove(0);
					for (j = 0; j < blockItemTotal; j++)
					{
						if (vectorArr[j] == 0)
						{
							if (linkage(dataArr.get(temp), dataArr.get(j)) <= threshold)
							{
								tempArr.add(j);
								vectorArr[j]	= n;
								edgeArr.add(new int[]{blockRowArr.get(temp), blockRowArr.get(j)});
							}
						}
					}
				}
			}
		}
		
		dataArr.clear();
		tempArr.clear();
		vectorArr	= null;
	}
	

	/*
	 * find clusters as connected components in a graph where edges are connection among records and vertices are records
	 */
	private ArrayList<Integer> findConnComp()
	{
		int i, rootU, rootV, edgeTotal, pointExactTotal;
		ArrayList<Integer> parentArr, weightArr;
	
		pointExactTotal	= clusterExactIndArr.size();
		parentArr		= new ArrayList<Integer>(pointExactTotal);
		weightArr		= new ArrayList<Integer>(pointExactTotal);
		
		for(i = 0; i < pointExactTotal; ++i)
		{
			parentArr.add(i);
			weightArr.add(0);
		}
	
		edgeTotal	= edgeArr.size();
	
		for(i = 0; i < edgeTotal; ++i)
		{
			rootU	= findRoot(edgeArr.get(i)[0], parentArr);
			rootV	= findRoot(edgeArr.get(i)[1], parentArr);
	
			if(rootU != rootV)
				makeEquivalent(rootU, rootV, parentArr, weightArr);
		}
	
		weightArr.clear();
		
		return parentArr;
	}
	
	
	/*
	 * find root of a point in components  
	 */
	private int findRoot(int pointID, ArrayList<Integer> parentArr)
	{
		if(parentArr.get(pointID) != pointID)
			parentArr.set(pointID, findRoot(parentArr.get(pointID), parentArr));
		return parentArr.get(pointID);
	}
	
	
	/*
	 * unify two components rooted at rootU and rootV
	 */
		
	private void makeEquivalent(int rootU, int rootV, ArrayList<Integer> parentArr, ArrayList<Integer> weightArr)
	{
		if(weightArr.get(rootU) < weightArr.get(rootV))
			parentArr.set(rootU, rootV);
		else if(weightArr.get(rootU) > weightArr.get(rootV))
			parentArr.set(rootV, rootU);
		else
		{
			parentArr.set(rootV, rootU);
			weightArr.set(rootU, weightArr.get(rootU) + 1);
		}
	}

	
	/*
	 * combine exact clusters and approximate clusters
	 */
	private void findFinalCluster(ArrayList<Integer> rootArr)
	{
		int i, j, indCluster;
		
		clusterArr	= new ArrayList<ArrayList<ArrayList<String>>>(rootArr.size());
		
		for(i = 0; i < rootArr.size(); ++i)
			clusterArr.add(new ArrayList<ArrayList<String>>());
		
		for(i = 0; i < rootArr.size(); ++i)
		{
			indCluster	= rootArr.get(i);

			for(j = 0; j < clusterExactIndArr.get(i).size(); ++j)
				clusterArr.get(indCluster).add(pointArr.get(clusterExactIndArr.get(i).get(j)));
		}
		
		clusterExactIndArr.clear();
		rootArr.clear();
	}
	
	
	/*
	 * get index of lastname attribute  
	 */
	private int getBlockInd()
	{
		int t = 0;
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		try 
		{
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(new File(configFileStr));
			Element blockEle	= (Element) doc.getElementsByTagName("block").item(0);
			Element blockIndEle	= (Element) blockEle.getElementsByTagName("index").item(0);
			String blockInd 	= blockIndEle.getChildNodes().item(0).getNodeValue().trim();
			t =  Integer.parseInt(blockInd);
			
//			Element blockTypeEle	= (Element) blockEle.getElementsByTagName("type").item(0);
//			String blockType	 	= blockTypeEle.getChildNodes().item(0).getNodeValue().trim();
//			alphabetType	=  Integer.parseInt(blockType);
			
			Element blockLenEle	= (Element) blockEle.getElementsByTagName("length").item(0);
			String blockLen 	= blockLenEle.getChildNodes().item(0).getNodeValue().trim();
			lmerSize	=  Integer.parseInt(blockLen) - lmerTotal + 1;
		} 
		catch (IOException e1) 
		{
			e1.printStackTrace();
		} 
		catch (SAXException e1) 
		{
			e1.printStackTrace();
		} 
		catch (ParserConfigurationException e1) 
		{
			e1.printStackTrace();
		}
		
		return t;
	}
	
	
	/*
	 * calculate linkage value (total distance) between two array of attributes
	 */
	private int linkage(ArrayList<String> a, ArrayList<String> b)
	{
		int w = 0;
		for (int i = 0; i < compArr.size(); i++)
		{
			if (compArr.get(i) == 1)
				w += EditDist(a, b, editAtt, threshold-w);
			else if (compArr.get(i) == 2)
				w += ReversalDist(a, b, reversalAtt, threshold-w);
			else if (compArr.get(i) == 3)
				w += TruncDist(a, b, truncAtt, threshold-w);			
		}
		
		return w;
	}

	
	/*
	 * edit distance between two array of attributes
	 */
	private int EditDist(ArrayList<String> a, ArrayList<String> b, ArrayList<int[]> compareAtt, int d)
	{
		int w, a_f, b_f, a_index, b_index, temp;
		int[] att;
		String s1, s2;
			
		w 	= 0;
		
		for (int g = 0; g < compareAtt.size(); g++)
		{
			att 	= compareAtt.get(g);
			a_f	= Integer.parseInt(a.get(a.size() - 1));
			b_f	= Integer.parseInt(b.get(b.size() - 1));
			a_index = att[a_f];
			b_index = att[b_f];
			s1 		= a.get(a_index);
			s2 		= b.get(b_index);
			temp 	= BasicED(s1,s2,d);
			w 		+= temp;
			d 		= d-temp;
		}
		
		return w;
	}
	
	
	/*
	 * calculate edit distance between two strings 
	 */
	private int BasicED(String str1, String str2, int d)
	{
		int dist = d;
		
		if (Math.abs(str1.length()-str2.length()) > dist)
			return (threshold + 1);
		else if (str1.equals(str2)) 
			return 0;
		else if (dist == 0)
			return (threshold + 1);
		else if (2*dist+1 >= Math.max(str1.length(), str2.length()))
			return BasicED2(str1,str2,dist);
		else
		{
			String s1, s2;
			if (str1.length() > str2.length())
			{
				s1 = str2;
				s2 = str1;
			}
			else
			{
				s1 = str1;
				s2 = str2;
			}
			int row = s1.length()+1;
			int col = 2*dist+1;
			int diagonal = dist + s2.length() - s1.length();
			int [][] m = new int[row][col];
			for (int i = 0; i < dist + 1; i++)
			{
				for (int j = dist-i; j < col; j++)
				{
					if (i == 0)
					{
						m[i][j] = j-dist;
					}
					else if (j == dist-i)
					{
						m[i][j] = m[i-1][j+1] + 1;
					}
					else if (j != col - 1)
					{
						if (s1.charAt(i-1) == s2.charAt(j-(dist-i)-1))
						{
							m[i][j] = Math.min(Math.min(m[i-1][j],m[i-1][j+1]+1),m[i][j-1]+1);
						}
						else
						{
							m[i][j] = Math.min(Math.min(m[i-1][j]+1,m[i-1][j+1]+1),m[i][j-1]+1);
						}
					}
					else 
					{
						if (s1.charAt(i-1) == s2.charAt(j-(dist-i)-1))
						{
							m[i][j] = Math.min(m[i-1][j],m[i][j-1]+1);
						}
						else
						{
							m[i][j] = Math.min(m[i-1][j]+1,m[i][j-1]+1);
						}
					}
					if (j == diagonal && m[i][j] > dist)
					{
						return (threshold + 1);
					}
				}
			}
			for (int i = dist + 1; i < s2.length() - dist + 1; i++)//(dist+1) + (s2.length() - (2*dist+1) + 1)
			{
				for (int j = 0; j < col; j++)
				{
					if (j == 0)
					{
						if (s1.charAt(i-1) == s2.charAt(j+(i-dist)-1))
						{
							m[i][j] = Math.min(m[i-1][j],m[i-1][j+1]+1);
						}
						else
						{
							m[i][j] = Math.min(m[i-1][j]+1,m[i-1][j+1]+1);
						}
					}
					else if (j != col - 1)
					{
						if (s1.charAt(i-1) == s2.charAt(j+(i-dist)-1))
						{
							m[i][j] = Math.min(Math.min(m[i-1][j],m[i-1][j+1]+1),m[i][j-1]+1);
						}
						else
						{
							m[i][j] = Math.min(Math.min(m[i-1][j]+1,m[i-1][j+1]+1),m[i][j-1]+1);
						}
					}
					else
					{
						if (s1.charAt(i-1) == s2.charAt(j+(i-dist)-1))
						{
							m[i][j] = Math.min(m[i-1][j],m[i][j-1]+1);
						}
						else
						{
							m[i][j] = Math.min(m[i-1][j]+1,m[i][j-1]+1);
						}
					}
					if (j == diagonal && m[i][j] > dist)
					{
						return (threshold + 1);
					}
				}
			}
			for (int i = s2.length() - dist + 1; i < row; i++)
			{
				for (int j = 0; j < col - i + s2.length() - dist; j++)//j<col - (i - (s2.length() - dist))
				{
					if (j == 0)
					{
						if (s1.charAt(i-1) == s2.charAt(j+(i-dist)-1))
						{
							m[i][j] = Math.min(m[i-1][j],m[i-1][j+1]+1);
						}
						else
						{
							m[i][j] = Math.min(m[i-1][j]+1,m[i-1][j+1]+1);
						}
					}
					else
					{
						if (s1.charAt(i-1) == s2.charAt(j+(i-dist)-1))
						{
							m[i][j] = Math.min(Math.min(m[i-1][j],m[i-1][j+1]+1),m[i][j-1]+1);
						}
						else
						{
							m[i][j] = Math.min(Math.min(m[i-1][j]+1,m[i-1][j+1]+1),m[i][j-1]+1);
						}
					}
					if (j == diagonal && m[i][j] > dist)
					{
						return (threshold + 1);
					}
				}
			}
			return (m[row-1][diagonal]);
		}
	}

	
	private int BasicED2(String s1, String s2, int d)
	{
		int row = s1.length()+1;
		int col = s2.length()+1;
		int [][] m = new int[row][col];
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
			{
				if (i == 0)
				{
					m[i][j] = j;
				}
				else if (j == 0)
				{
					m[i][j] = i;
				}
				else
				{
					if (s1.charAt(i-1) == s2.charAt(j-1))
					{
						m[i][j] = Math.min(Math.min(m[i-1][j-1],m[i-1][j]+1),m[i][j-1]+1);
					}
					else
					{
						m[i][j] = Math.min(Math.min(m[i-1][j-1]+1,m[i-1][j]+1),m[i][j-1]+1);
					}
				}
				if ((row-col) == (i-j) && (m[i][j] > d))
				{
					return (threshold + 1);
				}
			}
		}
		return (m[row-1][col-1]);
	}
	

	
	private int ReversalDist(ArrayList<String> a, ArrayList<String> b, ArrayList<ArrayList<int[]>> compareAtt, int d)
	{
		int w = 0;
		for (int g = 0; g < compareAtt.size(); g++)
		{
			ArrayList<int[]> att = compareAtt.get(g);
			int w1 = EditDist(a, b, att, d-w);
			ArrayList<int[]> attRev	= new ArrayList<int[]>();
			attRev.add(att.get(1));
			attRev.add(att.get(0));
			int w2 = EditDist(a ,b, attRev, d-w);
			w += Math.min(w1, w2);
		}
		return w;
	}

	private int TruncDist(ArrayList<String> a, ArrayList<String> b, ArrayList<int[]> compareAtt, int d)
	{
		int w = 0;
		ArrayList<String> a_ = new ArrayList<String>(a);
		ArrayList<String> b_ = new ArrayList<String>(b);

		int a_f = Integer.parseInt(a_.get(a_.size()-1));
		int b_f = Integer.parseInt(b_.get(b_.size()-1));
		
		for (int i = 0; i < compareAtt.size(); i++)
		{
			int ai = compareAtt.get(i)[a_f];
			String tempA = a_.get(ai);
			int bi = compareAtt.get(i)[b_f];
			String tempB = b_.get(bi);
			if (tempA.length() >= truncRate.get(i) && tempB.length() >= truncRate.get(i))
			{
				a_.set(ai, tempA.substring(0, truncRate.get(i)));
				b_.set(bi, tempB.substring(0, truncRate.get(i)));
			}
			else
			{
				int len = Math.min(tempA.length(), tempB.length());
				a_.set(ai, tempA.substring(0, len));
				b_.set(bi, tempB.substring(0, len));
			}
			
			//System.out.println(a_.get(ai) + "\t" + b_.get(bi));
		}
		
		
		w = EditDist(a_, b_, compareAtt, d);
		return w;
	}
		
	
	// read configuration file
	private String getTitle()
	{
		String fieldStr	= "";
		DocumentBuilderFactory dbf 		= DocumentBuilderFactory.newInstance();
		try 
		{
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(new File(configFileStr));
			Element elem = (Element)doc.getElementsByTagName("title").item(0);
			fieldStr	= elem.getChildNodes().item(0).getNodeValue().trim();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		return fieldStr;
	}
	private ArrayList<String> getInputFileNameList()
	{
		ArrayList<String> fileNameArr	= new ArrayList<String>();
		DocumentBuilderFactory dbf 		= DocumentBuilderFactory.newInstance();
		try {
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(new File(configFileStr));
			NodeList datasetList = doc.getElementsByTagName("dataset");
			for (int i = 0; i < datasetList.getLength(); i++)
			{
				Element dataset = (Element) datasetList.item(i);
				NodeList valuelist = dataset.getElementsByTagName("value");
				Element value = (Element) valuelist.item(0);
				NodeList childNodes = value.getChildNodes();
				for (int j = 0; j < childNodes.getLength(); j++)
				{
					String temp = childNodes.item(j).getNodeValue().trim();
					fileNameArr.add(temp);
				}
			}
		} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
		} catch (SAXException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();

		} catch (ParserConfigurationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		return fileNameArr;
	}
	
	private ArrayList<Integer> getCommonAttrList()
	{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		try {
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(new File(configFileStr));
			NodeList datasetList = doc.getElementsByTagName("dataset");
			for (int i = 0; i < datasetList.getLength(); i++)
			{
				ArrayList<Integer> tmp = new ArrayList<Integer>();
				Element index = (Element) datasetList.item(i);
				NodeList valuelist = index.getElementsByTagName("index");
				Element value = (Element) valuelist.item(0);
				NodeList childNodes = value.getChildNodes();
				String [] temp = childNodes.item(0).getNodeValue().trim().split(" ");
				for (int j = 0; j < temp.length; j++)
				{
					tmp.add(Integer.parseInt(temp[j]));
				}
				indexDatasetArr.add(tmp);
			}
		} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
		} catch (SAXException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();

		} catch (ParserConfigurationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		ArrayList<Integer> attrCommonArr	= new ArrayList<Integer>();
		
		for (int i = 0; i < indexDatasetArr.get(0).size(); ++i)
		{
			int j;
			for (j = 0; j < indexDatasetArr.size(); ++j)
			{
				if (indexDatasetArr.get(j).get(i) < 0)
					break;
			}
			
			if (j >= indexDatasetArr.size())
				attrCommonArr.add(i);
		}
		
		return attrCommonArr;
	}
	
	private int getInputThreshold()
	{
		int t = 0;
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		try {
			DocumentBuilder db 		= dbf.newDocumentBuilder();
			Document doc 			= db.parse(new File(configFileStr));
			Element thresholdVal 	= (Element) doc.getElementsByTagName("threshold").item(0);
			String threshold 		= thresholdVal.getChildNodes().item(0).getNodeValue().trim();
			t =  Integer.parseInt(threshold);
		} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
		} catch (SAXException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();

		} catch (ParserConfigurationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		return t;
	}
	
	public void getInputComparisonPara()
	{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		try {
			DocumentBuilder db = dbf.newDocumentBuilder();
			Document doc = db.parse(new File(configFileStr));
			NodeList comparisonList = doc.getElementsByTagName("comparison");
			for (int i = 0; i < comparisonList.getLength(); i++)
			{
				Element comparison = (Element) comparisonList.item(i);
				Element methodVal = (Element) comparison.getElementsByTagName("dist_calc_method").item(0);
				String comFun = methodVal.getChildNodes().item(0).getNodeValue().trim();
				if (comFun.equals("1") ||  comFun.equals("3"))//edit dist, truncate dist
				{
					Element comAttVal	= (Element) comparison.getElementsByTagName("comparing_attribute_indices").item(0);
					String att 			= comAttVal.getChildNodes().item(0).getNodeValue().trim();
					int [] attInt 		= new int[indexDatasetArr.size()];
					int j, attrInd;
					attrInd				= Integer.parseInt(att);
					//System.out.println("coomon: " + attrCommonArr);
					//System.out.println("indexDatasetArr: " + indexDatasetArr);
					for (j = 0; j < indexDatasetArr.size(); j++)
					{
						attInt[j] = indexDatasetArr.get(j).get(attrCommonArr.get(attrInd));
					}
					
					if (comFun.equals("1"))
					{
						addEditAtt(attInt);
					}
					else
					{
						Element rateVal = (Element) comparison.getElementsByTagName("truncation_count").item(0);
						String rate = rateVal.getChildNodes().item(0).getNodeValue().trim();
						int rateInt = Integer.parseInt(rate);
						addTruncAtt(attInt, rateInt);
					}
				}
				else if (comFun.equals("2"))
				{
					Element comAttVal1 = (Element) comparison.getElementsByTagName("comparing_attribute_indices").item(0);
					String [] att = comAttVal1.getChildNodes().item(0).getNodeValue().trim().split(",");					
					ArrayList<int[]> tempArr	= new ArrayList<int[]>();
					int j, k;
					for (k = 0; k < 2; ++k)
					{
						int [] attInt = new int[indexDatasetArr.size()];
						for (j = 0; j < indexDatasetArr.size(); j++)
							attInt[j] = indexDatasetArr.get(j).get(attrCommonArr.get(Integer.parseInt(att[k])));
							
						tempArr.add(attInt);
					}
					
					addReversalAtt(tempArr);
					
				}	
				else
				{
					//System.out.println("wrong chosen of comparison");
				}
			}
		} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
		} catch (SAXException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();

		} catch (ParserConfigurationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
	
	private void addEditAtt(int[] att)
	{
		if (!compArr.contains(1))
		{
			compArr.add(1);
		}
		editAtt.add(att);
	}

	private void addReversalAtt(ArrayList<int[]> att)
	{
		if (!compArr.contains(2))
		{
			compArr.add(2);
		}
		reversalAtt.add(att);
	}
	
	private void addTruncAtt(int[] att, int r)
	{
		if (!compArr.contains(3))
		{
			compArr.add(3);
		}
		truncAtt.add(att);
		truncRate.add(r);
	}
	
	// format and write output to external file
	
	private void output()
	{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		try 
		{
			DocumentBuilder db 	= dbf.newDocumentBuilder();
			Document doc 		= db.parse(new File(configFileStr));
			Element output 		= (Element) doc.getElementsByTagName("output_function").item(0);		
			Element outFileVal 	= (Element) output.getElementsByTagName("output_file").item(0);
			String outFile 		= outFileVal.getChildNodes().item(0).getNodeValue().trim();
			writeOutput(outFile);					
			
		} 
		catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
		} 
		catch (SAXException e1) 
		{
			// TODO Auto-generated catch block
			e1.printStackTrace();

		} 
		catch (ParserConfigurationException e1) 
		{
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
	
	private void writeOutput(String dirOutput)
	{
		try 
		{
			int i, j;
			BufferedWriter writer;
			ArrayList<String> recordRow;		
			int fileInd;
			
			String outFile		= dirOutput + "output.csv";
			//String outInfoFile	= dirOutput + "infoOutput.txt";
			
			clusterTotal	= 0;
			//System.out.println(outFile);	
			
			ArrayList<String> fileNameShortArr	= new ArrayList<String>();
			for (i = 0; i < fileNameArr.size(); ++i)
			{
				String[] nameSplitArr	= fileNameArr.get(i).split("\\/");
				fileNameShortArr.add(nameSplitArr[nameSplitArr.length - 1]);
			}
			
			writer			= new BufferedWriter(new FileWriter(outFile));
			writer.write("ClusterID,FileName," + fieldStr + "\n");
			for (i = 0; i < clusterArr.size(); i++)
			{
				if(clusterArr.get(i).size() <= 0)
					continue;
				++clusterTotal;								
				for(j = 0; j < clusterArr.get(i).size(); ++j)
				{
					recordRow	= clusterArr.get(i).get(j);
					fileInd		= Integer.parseInt(recordRow.get(recordRow.size() - 1));
					writer.write(clusterTotal + "," + fileNameShortArr.get(fileInd) + ",");
					for (int k = 0; k < indexDatasetArr.get(fileInd).size(); ++k)
						if (indexDatasetArr.get(fileInd).get(k) >= 0)
							writer.write(recordRow.get(indexDatasetArr.get(fileInd).get(k)) + ",");
						else
							writer.write(",");
					writer.write("\n");
				}		
			}
			writer.close();
			
			//writer			= new BufferedWriter(new FileWriter(outInfoFile));
			//writer.write(pointTotal + "\r\n");
			//writer.write(clusterTotal + "\r\n");
			//writer.write((System.currentTimeMillis() - startT) + "");
			//writer.close();
		} 
		catch (IOException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	public static boolean isNumeric(String str)  
	{  
	  try  
	  {  
	    double d = Double.parseDouble(str);  
	  }  
	  catch(NumberFormatException nfe)  
	  {  
	    return false;  
	  }  
	  return true;  
	}
}
