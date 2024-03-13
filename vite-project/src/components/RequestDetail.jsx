import { useState, useEffect } from "react";
import { useParams } from "react-router-dom";
import Table from "react-bootstrap/Table";
import React from "react";
import "./RequestDetail.css";

function RequestDetail() {
  let { id } = useParams();
  const [details, setDetails] = useState();
  const [result, setResult] = useState();
  const [detailId, setDetailId] = useState(id);

  useEffect(() => {
    fetch_details();
  }, [detailId]);

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

  function handleUpdateRequest(data) {
    setRequest(data);
    console.log(data);
  }

  function handleUpdateResult(data) {
    setResult(data);
    console.log(data);
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
              <td  className="dna">{details.alignment_request.dna_string}</td>
              <td>{details.alignment_request.status}</td>
              <td>{details.alignment_request.date_submitted}</td>
              <td>{details.alignment_request.date_updated}</td>
            </tr>
          </tbody>
        </Table>
      )}
      {details && details.alignment_results.hasOwnProperty("protein_id") && (
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
              <td className="dna">{details.alignment_results.protein_dna_seq}</td>
              <td>{details.alignment_results.organism}</td>
              <td>{details.alignment_results.filename}</td>
            </tr>
          </tbody>
        </Table>
      )}
    </>
  );
}

export default RequestDetail;
