import { useState, useEffect } from "react";
import { useParams } from "react-router-dom";
import Table from "react-bootstrap/Table";
import React from "react";
import { useInterval } from "useHooks-ts";
import "./RequestDetail.css";

function RequestDetail() {
  let { id } = useParams();
  const [details, setDetails] = useState();
  const [result, setResult] = useState();
  const [detailId, setDetailId] = useState(id);
  const [showWarning, setShowWarning] = useState(false);
  const [showResults, setShowResults] = useState(false);

  useEffect(() => {
    fetch_details();
  }, [detailId]);

  useEffect(() => {
    const sw =
      details &&
      details.alignment_request.status === "COMPLETE" &&
      !details.hasOwnProperty("alignment_results");
    setShowWarning(sw);
    const sr = details && details.hasOwnProperty("alignment_results");
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
        console.log(data);
        setDetails(data);
      });
  }

  function format_alignment_detail(details) {
    if (details) {
      lines = details.split("\n");
      return (
        <>
          {lines.map((line) => {
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
                <td>{details.alignment_results.protein_id}</td>
                <td className="details">
                  {format_alignment_detail(
                    details.alignment_results.alignment_detail
                  )}
                </td>
                <td className="table-dna">
                  {details.alignment_results.protein_dna_seq}
                </td>
                <td>{details.alignment_results.organism}</td>
                <td>{details.alignment_results.filename}</td>
              </tr>
            </tbody>
          </Table>
        </div>
      )}
    </>
  );
}

export default RequestDetail;
