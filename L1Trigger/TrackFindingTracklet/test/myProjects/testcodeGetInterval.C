


double getIntervalContainingFractionOfEntries(TH1* histogram, double interval, int minEntries = 5);

void testcodeGetInterval () {
    TH1F* h = new TH1F("h", "h", 10, -5, 5);

    h->Fill(1);
    h->Fill(1);
    h->Fill(1);
    h->Fill(3);
    h->Fill(3);
    h->Fill(3);
    h->Fill(20);
    h->Fill(20);
    h->Fill(20);
    h->Fill(20);
    
    
    

    // h->Fill(-4.9);
    // h->Fill(-4);
    // h->Fill(-3);
    // h->Fill(-2);
    // h->Fill(-1);
    // h->Fill(0);
    // h->Fill(1);
    // h->Fill(2);
    // h->Fill(3);
    // h->Fill(4.5);



    double interval = getIntervalContainingFractionOfEntries(h, 0.6);

    cout << "interval: " << interval << endl;
    cout << "overflow bin entries: " <<  h->GetBinContent(21) << endl;

    h->Draw();
}


double getIntervalContainingFractionOfEntries(TH1* absResidualHistogram, double quantileToCalculate, int minEntries) {
  double totalIntegral = absResidualHistogram->Integral(0, absResidualHistogram->GetNbinsX() + 1);
  double numEntries = absResidualHistogram->GetEntries();

  // Check that the interval is not somewhere in the overflow bin
  double maxAllowedEntriesInOverflow = totalIntegral * (1 - quantileToCalculate);
  double nEntriesInOverflow = absResidualHistogram->GetBinContent(absResidualHistogram->GetNbinsX() + 1);
  if (nEntriesInOverflow > maxAllowedEntriesInOverflow) {
    // cout << "WARNING : Cannot compute range corresponding to interval, as it is in the overflow bin" << endl;
    return absResidualHistogram->GetXaxis()->GetXmax() * 1.2;
  }

  // Calculate quantile for given interval
  double interval[1];
  double quantile[1] = {quantileToCalculate};
  if (totalIntegral > 0 && numEntries >= minEntries) {
    absResidualHistogram->GetQuantiles(1, interval, quantile);
  } else {
    cout << "WARNING: histo " << absResidualHistogram->GetName()
         << " empty or with too few entries, so can't calc quantiles." << endl;
    interval[0] = 0.;
  }

  return interval[0];
}