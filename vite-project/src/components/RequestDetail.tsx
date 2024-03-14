import { useState, useEffect } from "react";
import { useParams } from "react-router-dom";
import Table from "react-bootstrap/Table";
import { useInterval } from "useHooks-ts";
import "./RequestDetail.css";

type DetailsType = {
  alignment_request: AlignmentRequestType;
  alignment_results?: AlignmentResultType;
};

type AlignmentRequestType = {
  dna_string: string;
  status: string;
  date_submitted: string;
  date_updated: string;
};

type AlignmentResultType = {
  protein_id: string;
  alignment_detail: string;
  protein_dna_seq: string;
  organism: string;
  filename: string;
};

function RequestDetail() {
  let { id } = useParams();
  const [details, setDetails] = useState<DetailsType | undefined>();
  const [detailId, setDetailId] = useState(id);
  const [showWarning, setShowWarning] = useState<boolean | undefined>(false);
  const [showResults, setShowResults] = useState<boolean | undefined>(false);

  useEffect(() => {
    fetch_details();
  }, [detailId]);

  useEffect(() => {
    const sw: boolean | undefined =
      details &&
      details.alignment_request.status === "COMPLETE" &&
      !details.hasOwnProperty("alignment_results");
    setShowWarning(sw);
    const sr: boolean | undefined =
      details && details.hasOwnProperty("alignment_results");
    setShowResults(sr);
  }, [details]);

  // grabs the details every 30 seconds until the results load
  useInterval(
    () => {
      fetch_details();
    },
    details && details.hasOwnProperty("alignment_results") ? null : 30000
  );

  function fetch_details() {
    fetch("/api/detail/" + detailId, {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setDetails(data);
      });
  }

  function format_alignment_detail(details:string|undefined) {
    if (details) {
      const lines:Array<string> = details.split("\n");
      return (
        <>
          {lines.map((line:string) => {
            return <p className="line">{line}</p>;
          })}
        </>
      );
    }
  }

  return (
    <>
      {details && (
        <div className="table-squish">
          <Table striped bordered hover>
            <thead>
              <tr>
                <td>DNA Query</td>
                <td>Status</td>
                <td>Date Submitted</td>
                <td>Date Updated</td>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td className="table-dna">
                  {details.alignment_request.dna_string}
                </td>
                <td className="table-status">
                  {details.alignment_request.status}
                </td>
                <td className="table-date">
                  {details.alignment_request.date_submitted}
                </td>
                <td className="table-date">
                  {details.alignment_request.date_updated}
                </td>
              </tr>
            </tbody>
          </Table>
        </div>
      )}
      {showWarning && (
        <h3 className="no-match">No match was found for this query.</h3>
      )}
      {showResults && (
        <div className="table-squish">
          <Table striped bordered hover>
            <thead>
              <tr>
                <td>Protein id</td>
                <td>Alignment Result</td>
                <td>Protein DNA Target</td>
                <td>Organism</td>
                <td>Filename</td>
              </tr>
            </thead>
            <tbody>
              <tr>
                <td>{details?.alignment_results?.protein_id}</td>
                <td className="details">
                  {format_alignment_detail(
                    details?.alignment_results?.alignment_detail
                  )}
                </td>
                <td className="table-dna">
                  {details?.alignment_results?.protein_dna_seq}
                </td>
                <td>{details?.alignment_results?.organism}</td>
                <td>{details?.alignment_results?.filename}</td>
              </tr>
            </tbody>
          </Table>
        </div>
      )}
    </>
  );
}

export default RequestDetail;
