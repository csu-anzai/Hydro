


	/*====================================================================
		
		Create Gain Region Mask 
		This flags whether a cell is inside the gain region or 
		outside the gain region.

		The gain region is defined as the region bounded by:
		(start)	The point at which the lateral average of qenr 
				transitions from negative to positive
		(end)	The maximum of the shock radius and the point at 
				which the lateral average of qenr goes negative

		(1)	Average qenr
		(2)	Compute the radial positions where qenr goes from 
			(a)	negative to positive (qenrNP)
			(a)	positive to negative (qenrPN)
		(3)	For every ray
			(1)	For every cell, determine if the cell lies inside the
				gain region

				If it does, set the mask to 1.0
				Otherwise, set the mask to 0.0

	=======================================================================*/

	double qenrNP = 0.0;
	double qenrPN = 0.0;
	bool foundNP = false, foundPN = false;

	double gainRadius = 0.0;

	// Average qenr
	std::vector<std::string> shellVars;
	shellVars.push_back("qenr");

	vtkSmartPointer<vtkDataArray> radCoords = vtkRectilinearGrid::SafeDownCast(data.getDataSet())->GetXCoordinates();
	ShellAveragePlaneOperator shellAvg(	radCoords,
										XAXIS,
										CS_SPHERE,
										shellVars	);
	shellAvg.process(data.getDataSet());

	vtkSmartPointer<vtkDataArray> qenrAvg = shellAvg.getDataArray("qenr");

	// Find qenrNP
	// For every bin in the interior (not on the edge) of the domain,
	// look at the previous and next bin values.
	// If the next is positive and the previous negative, interpolate
	// to find the coordinate where qenr becomes zero
	double yPrev = 0.0;
	double yNext = 1.0;
	for(int i=qenrAvg->GetNumberOfTuples()-1; i>1; --i)
	{
		yPrev = qenrAvg->GetTuple1(i);
		yNext = qenrAvg->GetTuple1(i-1);


		if( ((yPrev < 0.0) && (yNext > 0.0)) || ((yPrev > 0.0) && (yNext < 0.0)) )
		{
			// Get the volumetric centers of the previous and next bins
			double xPrevBndL = radCoords->GetTuple1(i-1);
			double xPrevBndR = radCoords->GetTuple1(i);
			double xNextBndL = radCoords->GetTuple1(i-2);
			double xNextBndR = radCoords->GetTuple1(i-1);
			double xPrev = pow(0.5*(xPrevBndL*xPrevBndL*xPrevBndL + xPrevBndR*xPrevBndR*xPrevBndR), 1.0/3.0);
			double xNext = pow(0.5*(xNextBndL*xNextBndL*xNextBndL + xNextBndR*xNextBndR*xNextBndR), 1.0/3.0);
			double slope = (yNext-yPrev)/(xNext-xPrev);
			double xCrossing = -yPrev/slope + xPrev;

			// // qenr is going from negative to positive => qenrNP = xCrossing
			// if( (yPrev < 0.0) && (yNext > 0.0) && !foundNP)
			// {
			// 	qenrNP = xCrossing;
			// 	foundNP = true;
			// }
			// // qenr is going from positive to negative => qenrPN = xCrossing
			// else if( (yPrev > 0.0) && (yNext < 0.0) && !foundPN )
			// {
			// 	qenrPN = xCrossing;
			// 	foundPN = true;
			// }

			if( (yPrev > 0.0) && (yNext < 0.0) )
			{
				qenrPN = xCrossing;
				break;
			}
		}
	}

	// Set the gain radius as the transition from negative to positive
	gainRadius = qenrPN; 

	// Store
	dataList.addScalar(gainRadius,"gainradius");

	// Compute the MFR through the mean gain radius surface
	double MFRgainRadius = 0.0;
	for (int i=0; i<data.getDataSet()->GetNumberOfCells(); i++)
	{
		// Determine the area at the volume centroid of the cell that the radial flow sees
		double cellEdges[6];
		data.getDataSet()->GetCellBounds(i,cellEdges);

		// std::cout<<cellEdges[0]<<"\t"<<cellEdges[1]<<"\t"<<gainRadius<<std::endl;

		// Determine if the cell is cut by the gain radius
		if(cellEdges[0]<gainRadius && cellEdges[1]>=gainRadius)
		{
			// std::cout<<"In here!"<<std::endl;

			double cellArea = 1.0;
			cellArea *= gainRadius*gainRadius; // Fixed radius
			cellArea *= -cos(cellEdges[3]) + cos(cellEdges[2]); // Phi contribution
			if(nDim>2)
			{
				cellArea *= cellEdges[5]-cellEdges[4]; // theta contribution
			}
			else
			{
				cellArea *= 2.0*vtkMath::Pi();
			}

			// Correct for the full 4pir^2 area

			cellArea *= volumeCorrection;

			MFRgainRadius += data["dens"]->GetTuple1(i)*data["velx"]->GetTuple1(i)*cellArea;
		}
	}

	// Store
	dataList.addScalar(MFRgainRadius,"mfrGR");

	std::cout<<"Gain radius: "<<gainRadius<<std::endl;
	std::cout<<"MFR at gain radius: "<<MFRgainRadius/MSOLAR<<std::endl;



	// Construct mask
	data.addVariable("mask");
	vtkSmartPointer<vtkDataArray> mask = data["mask"];

	for(int p = 0; p<meshDimensions[1]; ++p)
	{
		for(int t = 0; t< (nDim>2 ? meshDimensions[2] : 1); ++t)
		{
			startingIndex[XAXIS]	= 0;
			startingIndex[YAXIS]	= p;
			startingIndex[ZAXIS]	= t;
			std::vector<int> cellIndices = data.getCellsAlongRay(XAXIS, startingIndex);

			bool foundFirst = false, foundLast = false;
			double shockPos = 0.0;

			for (std::vector<int>::reverse_iterator i = cellIndices.rbegin(); i != cellIndices.rend(); ++i)
			{
				
				bool isShocked = shock->GetTuple1(*i) > 0.5; // The data is either 0 or 1, so you need to pick something reasonable

				// If it's shocked, add the radial coordinate to the sum
				if(isShocked)
				{
					if(!foundFirst)
						foundFirst = true;
				}
				else
				{
					if(foundFirst && !foundLast)
					{
						foundLast = true;
						shockPos = radius->GetTuple1(*i);
						// std::cout<<"shock pos = "<<shockPos<<std::endl;
					}
				}

				// Get cell coordinates
				double cellEdges[6];
				data.getDataSet()->GetCellBounds(*i,cellEdges);

				double xMin = cellEdges[0];
				double xMax = cellEdges[1];
				double maskValue = 0.0;

				// if(xMax > qenrNP)
				// {
				// 	if(xMin < std::max(qenrPN, shockPos))
				// 	{
				// 		maskValue = 1.0;
				// 	}
				// 	else
				// 	{
				// 		maskValue = 0.0;
				// 	}
				// }
				// else
				// {
				// 	maskValue = 0.0;
				// }

				if(	xMin>gainRadius && 
				   	xMin<shockPos &&
				   	xMax>gainRadius &&
				   	xMax<shockPos)
				{
					maskValue = 1.0;
				}

				mask->SetTuple1(*i,maskValue);
			}
		}
	}

	// Create mass variable
	std::clog<<"Creating mass variable..."; std::clog.flush();
	data.storeData(data["dens"]*data["cellVolume"],"mass");
	double Mgain = data.sum("mass","mask");

	std::clog<<"done!"<<std::endl;

	std::cout<<"Total species mass [Msun]: "<<data.sum("mass")/MSOLAR<<std::endl;
	std::cout<<"Gain region mass:"<<std::endl;
	std::cout<<"\t    Mass [g]: "<<Mgain<<std::endl;
	std::cout<<"\t Mass [Msun]: "<<Mgain/MSOLAR<<std::endl;


	// Store
	dataList.addScalar(Mgain,"massgainregion");
