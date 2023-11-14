


void testcodeQuantile () {
    TH1F* h = new TH1F("h", "h", 20, -10, 10);

    h->Fill(1);
    h->Fill(3);
    h->Fill(20);
    h->Fill(20);
    h->Fill(20);
    h->Fill(20);
    h->Fill(20);
    h->Fill(20);

    double interval[1];
    double quantile[1] = {0.75};

    h->GetQuantiles(1, interval, quantile);

    cout << "interval: " << interval[0] << endl;
    cout << "overflow bin entries: " <<  h->GetBinContent(21) << endl;

    h->Draw();
}